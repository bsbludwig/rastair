/// Run up to 5 closures/expressions in parallel via `rayon::join` and collect
/// results into a flat tuple.
#[macro_export]
macro_rules! rayon_all {
    ($a:expr $(,)?) => {
        $a
    };
    ($a:expr, $b:expr $(,)?) => {
        ::rayon::join(|| $a, || $b)
    };
    ($a:expr, $b:expr, $c:expr $(,)?) => {{
        let (a, (b, c)) = ::rayon::join(|| $a, || ::rayon::join(|| $b, || $c));
        (a, b, c)
    }};
    ($a:expr, $b:expr, $c:expr, $d:expr $(,)?) => {{
        let (a, (b, (c, d))) =
            ::rayon::join(|| $a, || ::rayon::join(|| $b, || ::rayon::join(|| $c, || $d)));
        (a, b, c, d)
    }};
    ($a:expr, $b:expr, $c:expr, $d:expr, $e:expr $(,)?) => {{
        let (a, (b, (c, (d, e)))) = ::rayon::join(
            || $a,
            || ::rayon::join(|| $b, || ::rayon::join(|| $c, || ::rayon::join(|| $d, || $e))),
        );
        (a, b, c, d, e)
    }};
}

#[cfg(test)]
mod tests {
    #[test]
    #[allow(unused, reason = "just testing the macro")]
    fn call() {
        fn id(x: i32) -> i32 {
            x
        }

        let one = rayon_all!(id(1));
        let (one, two) = rayon_all!(id(1), id(2));
        let (one, two, three) = rayon_all!(id(1), id(2), id(3));
        let (one, two, three, four) = rayon_all!(id(1), id(2), id(3), id(4));
        let (one, two, three, four, five) = rayon_all!(id(1), id(2), id(3), id(4), id(5));

        assert_eq!((one, two, three, four, five), (1, 2, 3, 4, 5));
    }
}
