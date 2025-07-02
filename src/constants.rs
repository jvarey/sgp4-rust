use std::f64::consts::PI;

// --- mathematical operations ---
pub const SMALL: f64 = 1e-10;

// distances
pub const KM2M: f64 = 1e3;
pub const FT2M: f64 = 0.3048;
pub const MILE2M: f64 = 1609.344;
pub const NM2M: f64 = 1852.0;
pub const MILE2FT: f64 = 5280.0;
pub const MILEPH2KMPH: f64 = 0.44704;
pub const NMPH2KMPH: f64 = 0.5144444;

// time
pub const DAY2SEC: f64 = 86400.0;
pub const DAY2MIN: f64 = 1440.0;
pub const DAY2HR: f64 = 24.0;
pub const HR2SEC: f64 = 3600.0;
pub const MIN2SEC: f64 = 60.0;
pub const YR2DAY: f64 = 365.25;
pub const CENT2YR: f64 = 100.0;
pub const CENT2DAY: f64 = CENT2YR * YR2DAY;

// angles
pub const HALFPI: f64 = PI / 2.0;
pub const TWOPI: f64 = 2.0 * PI;
pub const DEG2MIN: f64 = 60.0;
pub const DEG2ARCSEC: f64 = DEG2MIN * MIN2SEC;
pub const ARCSEC2RAD: f64 = 1.0 / DEG2ARCSEC * PI / 180.0;
pub const DEG2SEC: f64 = 360.0 / DAY2SEC;
pub const DEG2HR: f64 = 360.0 / DAY2HR;
pub const HR2RAD: f64 = DEG2HR * 1.0 * PI / 180.0;

// --- astrodynamics operations ---

// time
/// Julian date of the epoch J2000.0 (noon)
pub const J2000: f64 = 2451545.0;
/// Julian date of the epoch J2000.0 in UTC (midnight)
pub const J2000_UTC: f64 = 2451544.5;
/// offset between Julian dates and Modified Julian dates
pub const JD_TO_MJD_OFFSET: f64 = 2400000.5;

// EGM-08 (Earth) constants used here
pub const RE: f64 = 6378.1363;
pub const FLAT: f64 = 1.0 / 298.257223563;
pub const EARTHROT: f64 = 7.292115e-5;
pub const MU: f64 = 398600.4415;
pub const MUM: f64 = 3986004415e14;
pub const J2: f64 = 0.001082626174;
pub const J4: f64 = -1.6198976e-06;
