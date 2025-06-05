use crate::utils::Base;

#[derive(Debug, Default)]
pub struct Counter {
    pub c: usize,
    pub t: usize,
    pub a: usize,
    pub g: usize,
}

impl Counter {
    /// Interesting if there are multiple different bases seen
    pub fn multiple_bases(&self) -> bool {
        let mut count = 0;
        if self.c > 0 {
            count += 1;
        }
        if self.t > 0 {
            count += 1;
        }
        if self.a > 0 {
            count += 1;
        }
        if self.g > 0 {
            count += 1;
        }
        count >= 1
    }
}

impl FromIterator<Base> for Counter {
    fn from_iter<I: IntoIterator<Item = Base>>(iter: I) -> Self {
        let mut counter = Counter { c: 0, t: 0, a: 0, g: 0 };
        for c in iter {
            match c {
                Base::C => counter.c += 1,
                Base::T => counter.t += 1,
                Base::A => counter.a += 1,
                Base::G => counter.g += 1,
            }
        }
        counter
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counter_sanity() {
        let bases = vec![Base::A, Base::C, Base::G, Base::T, Base::A, Base::C];
        let counter: Counter = bases.into_iter().collect();
        assert_eq!(counter.a, 2);
        assert_eq!(counter.c, 2);
        assert_eq!(counter.g, 1);
        assert_eq!(counter.t, 1);
    }
}
