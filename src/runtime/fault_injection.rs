//! Injecting panics, errors and signals at fixed places, to test how the
//! pipelines handle them.
//!
//! Each environment variable takes a comma-separated list of `<point>` or
//! `<point>@<n>` (the first or n-th time, 1-based, that the named
//! [`FaultPoint`] is reached); signals are written as `<signal>:<point>[@<n>]`
//! with `<signal>` being `INT` or `TERM`. Only read in debug builds.

use color_eyre::eyre::Result;

/// Places where a fault can be injected, see [`fault_point`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaultPoint {
    /// A worker, before processing a segment
    Worker,
    /// The writer, after writing a segment
    Writer,
    /// The GPU inference thread, after a dispatch returned but before its
    /// results are handed out, i.e. while workers are waiting and queueing
    GpuDispatch,
}

impl FaultPoint {
    #[cfg(debug_assertions)]
    const ALL: [Self; 3] = [Self::Worker, Self::Writer, Self::GpuDispatch];

    /// Name used in the environment variables
    pub const fn name(self) -> &'static str {
        match self {
            Self::Worker => "worker",
            Self::Writer => "writer",
            Self::GpuDispatch => "gpu-dispatch",
        }
    }
}

/// Panic at the given points.
pub const INJECT_PANIC_VAR: &str = "RASTAIR_INJECT_PANIC";
/// Return an error from [`fault_point`] at the given points.
pub const INJECT_ERROR_VAR: &str = "RASTAIR_INJECT_ERROR";
/// Send a signal to the process (from the thread reaching the point), which
/// kills it.
pub const INJECT_SIGNAL_VAR: &str = "RASTAIR_INJECT_SIGNAL";

/// Inject whatever the environment asks for at this point. Compiled out of
/// release builds, where it always returns `Ok`.
#[inline]
#[track_caller]
pub fn fault_point(point: FaultPoint) -> Result<()> {
    #[cfg(debug_assertions)]
    return injection::reached(point);
    #[cfg(not(debug_assertions))]
    {
        let _ = point;
        Ok(())
    }
}

#[cfg(debug_assertions)]
mod injection {
    use super::{FaultPoint, INJECT_ERROR_VAR, INJECT_PANIC_VAR, INJECT_SIGNAL_VAR};
    use color_eyre::eyre::{Result, bail};
    use std::sync::{
        OnceLock,
        atomic::{AtomicUsize, Ordering},
    };
    use tracing::warn;

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) enum Fault {
        Panic,
        Error,
        Signal(i32),
    }

    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct Target {
        pub(super) fault: Fault,
        pub(super) point: FaultPoint,
        pub(super) nth: usize,
    }

    #[track_caller]
    #[expect(clippy::panic, reason = "injecting panics is the point")]
    pub(super) fn reached(point: FaultPoint) -> Result<()> {
        static TARGETS: OnceLock<Vec<(Target, AtomicUsize)>> = OnceLock::new();

        let targets = TARGETS.get_or_init(|| {
            [(INJECT_PANIC_VAR, Fault::Panic), (INJECT_ERROR_VAR, Fault::Error)]
                .into_iter()
                .map(|(var, fault)| (var, Some(fault)))
                .chain([(INJECT_SIGNAL_VAR, None)])
                .filter_map(|(var, fault)| Some((var, fault, std::env::var(var).ok()?)))
                .flat_map(|(var, fault, spec)| {
                    parse(&spec, fault)
                        .inspect_err(|error| warn!(%error, "Ignoring invalid {var}"))
                        .unwrap_or_default()
                })
                .map(|target| (target, AtomicUsize::new(0)))
                .collect()
        });

        for (target, reached) in targets.iter().filter(|(target, _)| target.point == point) {
            if reached.fetch_add(1, Ordering::Relaxed) + 1 != target.nth {
                continue;
            }
            match target.fault {
                Fault::Panic => panic!(
                    "Injected panic at {} (time {}) via {INJECT_PANIC_VAR}",
                    point.name(),
                    target.nth
                ),
                Fault::Error => bail!(
                    "Injected error at {} (time {}) via {INJECT_ERROR_VAR}",
                    point.name(),
                    target.nth
                ),
                Fault::Signal(signal) => raise(signal)?,
            }
        }
        Ok(())
    }

    #[cfg(unix)]
    fn raise(signal: i32) -> Result<()> {
        // SAFETY: raise(3) has no memory-safety preconditions.
        if unsafe { libc::raise(signal) } != 0 {
            bail!("Failed to raise injected signal: {}", std::io::Error::last_os_error());
        }
        Ok(())
    }

    #[cfg(not(unix))]
    fn raise(_signal: i32) -> Result<()> {
        bail!("Injecting signals is only supported on Unix")
    }

    /// `fault` is `None` for signal specs, which name their signal.
    pub(super) fn parse(spec: &str, fault: Option<Fault>) -> Result<Vec<Target>, String> {
        spec.split(',').map(|entry| parse_one(entry.trim(), fault)).collect()
    }

    fn parse_one(entry: &str, fault: Option<Fault>) -> Result<Target, String> {
        let (fault, entry) = match fault {
            Some(fault) => (fault, entry),
            None => {
                let (signal, entry) = entry
                    .split_once(':')
                    .ok_or_else(|| format!("expected `<signal>:<point>`, got `{entry}`"))?;
                (Fault::Signal(parse_signal(signal)?), entry)
            }
        };
        let (name, nth) = match entry.split_once('@') {
            Some((name, nth)) => {
                let nth = nth.parse().map_err(|error| format!("invalid count `{nth}`: {error}"))?;
                (name, nth)
            }
            None => (entry, 1),
        };
        if nth == 0 {
            return Err("the count is 1-based".to_string());
        }
        let point = FaultPoint::ALL.into_iter().find(|p| p.name() == name).ok_or_else(|| {
            let known = FaultPoint::ALL.map(FaultPoint::name).join(", ");
            format!("unknown fault point `{name}`, expected one of: {known}")
        })?;
        Ok(Target { fault, point, nth })
    }

    fn parse_signal(signal: &str) -> Result<i32, String> {
        match signal {
            #[cfg(unix)]
            "INT" => Ok(libc::SIGINT),
            #[cfg(unix)]
            "TERM" => Ok(libc::SIGTERM),
            _ => Err(format!("unsupported signal `{signal}`, expected INT or TERM")),
        }
    }
}

#[cfg(all(test, debug_assertions))]
mod tests {
    use super::FaultPoint;
    use super::injection::{Fault, Target, parse};

    #[test]
    fn parses_fault_specs() {
        let panic = Some(Fault::Panic);
        assert_eq!(
            parse("worker", panic),
            Ok(vec![Target { fault: Fault::Panic, point: FaultPoint::Worker, nth: 1 }])
        );
        assert_eq!(
            parse("gpu-dispatch@3", Some(Fault::Error)),
            Ok(vec![Target { fault: Fault::Error, point: FaultPoint::GpuDispatch, nth: 3 }])
        );
        assert!(parse("writer@0", panic).is_err());
        assert!(parse("writer@x", panic).is_err());
        assert!(parse("nope", panic).is_err_and(|error| error.contains("gpu-dispatch")));
    }

    #[cfg(unix)]
    #[test]
    fn parses_signal_specs() {
        assert_eq!(
            parse("TERM:worker@3, INT:writer", None),
            Ok(vec![
                Target { fault: Fault::Signal(libc::SIGTERM), point: FaultPoint::Worker, nth: 3 },
                Target { fault: Fault::Signal(libc::SIGINT), point: FaultPoint::Writer, nth: 1 },
            ])
        );
        assert!(parse("worker", None).is_err(), "signal specs name their signal");
        assert!(parse("KILL:worker", None).is_err());
    }
}
