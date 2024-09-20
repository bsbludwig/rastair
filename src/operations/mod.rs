use std::str::FromStr;

pub mod count_reads;
pub mod count_variants;

/// Represent a read softmask, to exclude certain portions of the read
#[derive(Clone, Copy, Debug)]
pub struct ReadMask(usize, usize);

#[derive(Clone, Copy, Debug)]
pub struct ReadMaskSetting
{
    pub r1: ReadMask,
    pub r2: ReadMask,
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseMaskError;

impl FromStr for ReadMaskSetting
{
    type Err = ParseMaskError;

    fn from_str(s: &str) -> Result<Self, Self::Err>
    {
        let values = s.trim().split(',').collect::<Vec<&str>>();
        if values.len() != 4 {
            return Err(ParseMaskError);
        }

        let r1_left = values[0].parse::<usize>().map_err(|_| ParseMaskError)?;
        let r1_right = values[1].parse::<usize>().map_err(|_| ParseMaskError)?;
        let r2_left = values[2].parse::<usize>().map_err(|_| ParseMaskError)?;
        let r2_right = values[3].parse::<usize>().map_err(|_| ParseMaskError)?;

        Ok(ReadMaskSetting { r1: ReadMask(r1_left, r1_right),
                             r2: ReadMask(r2_left, r2_right) })
    }
}

/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests
{

    // For testing
    use super::*;
    use anyhow::Result;
    use std::str::FromStr;

    #[test]
    fn can_create_from_string() -> Result<()>
    {
        let input_string = "1,2,3,4";
        let test = ReadMaskSetting::from_str(input_string).unwrap();

        assert_eq!(test.r1.0, 1);
        assert_eq!(test.r1.1, 2);
        assert_eq!(test.r2.0, 3);
        assert_eq!(test.r2.1, 4);
        Ok(())
    }

    #[test]
    fn tolerates_zero() -> Result<()>
    {
        let input_string = "1,0,0,1";
        let test = ReadMaskSetting::from_str(input_string).unwrap();

        assert_eq!(test.r1.0, 1);
        assert_eq!(test.r1.1, 0);
        assert_eq!(test.r2.0, 0);
        assert_eq!(test.r2.1, 1);
        Ok(())
    }

    #[test]
    fn fails_with_letters() -> Result<()>
    {
        let input_string = "1,b,a,1";
        let err = ReadMaskSetting::from_str(input_string).unwrap_err();
        assert_eq!(err, ParseMaskError);
        Ok(())
    }
}
