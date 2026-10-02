//! How a command runs as a process: its threads, how it handles a crash, how
//! segments reach the writer, how its output files appear, and how it reports
//! progress. Shared by `call`, `per-read` and `bam`.

pub mod fault_injection;
pub mod partial_output;
pub(crate) mod progress;
pub mod segments;
pub mod threads;
