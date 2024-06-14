// Bundle some utility structs and constants in here for all operations
pub const MAX_DEPTH: u32 = 10_000;

/// Utility structure to encode sam flags by name
#[allow(dead_code)]
pub struct Flags
{
    pub is_paired: u16,
    pub is_properly_paired: u16,
    pub is_unmapped: u16,
    pub mate_is_unmapped: u16,
    pub is_reverse_strand: u16,
    pub mate_is_reverse_strand: u16,
    pub is_first_in_pair: u16,
    pub is_second_in_pair: u16,
    pub is_not_primary: u16,
    pub is_failed: u16,
    pub is_duplicate: u16,
    pub is_supplemental: u16,
}
pub const FLAGS: Flags = Flags
{
    is_paired: 0x1,
    is_properly_paired: 0x2,
    is_unmapped: 0x4,
    mate_is_unmapped: 0x8,
    is_reverse_strand: 0x10,
    mate_is_reverse_strand: 0x20,
    is_first_in_pair: 0x40,
    is_second_in_pair: 0x80,
    is_not_primary: 0x100,
    is_failed: 0x200,
    is_duplicate: 0x400,
    is_supplemental: 0x800
};