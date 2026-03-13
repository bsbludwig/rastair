use crate::utils::Base;

#[derive(Debug, Default)]
pub struct Counter {
    pub c: usize,
    pub t: usize,
    pub a: usize,
    pub g: usize,
}

impl Counter {
    pub fn entries(&self) -> [(Base, usize); 4] {
        [(Base::C, self.c), (Base::T, self.t), (Base::A, self.a), (Base::G, self.g)]
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
                Base::Unknown => continue, // Ignore unknown bases
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
