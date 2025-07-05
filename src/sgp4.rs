use std::error::Error;

use crate::constants::*;
use crate::utils::*;

pub enum TypeRun {
    Catalog,
    Verification,
    FromJD,
    Manual,
}

#[derive(Default)]
pub enum Classification {
    #[default]
    Unclassified,
    Classified,
}

impl Classification {
    fn from(c: char) -> Option<Classification> {
        match c {
            'U' => Some(Classification::Unclassified),
            'C' => Some(Classification::Classified),
            _ => None,
        }
    }
}

pub enum PropagationError {
    InvalidElements,
    NegativeMeanMotion,
    EccentricityOutOfRange,
    NegativeSemilatusRectum,
    OrbitalDecay,
}

#[derive(Debug)]
pub struct TLE {
    line1: String,
    line2: String,
}

impl TLE {
    pub fn new(line1: String, line2: String) -> TLE {
        TLE { line1, line2 }
    }
}

#[derive(Default)]
struct SGP4InitOutput {
    ainv: f64,
    ao: f64,
    con41: f64,
    con42: f64,
    cosio: f64,
    cosio2: f64,
    eccsq: f64,
    omeosq: f64,
    posq: f64,
    rp: f64,
    rteosq: f64,
    sinio: f64,
    gsto: f64,
    no_unkozai: f64,
}

pub struct SGP4 {
    model: WGSModel,
    use_afspc_mode: bool,
    grav_const: GravitationalConstants,
    satrec: SatRec,
    use_deep_space: bool,
    x2ox3: f64,
    jdstart_full: f64,
    jdstop_full: f64,
    sgp4init_out: SGP4InitOutput,
    // TODO: implement deep space
    // ds: Option<...>
}

impl Default for SGP4 {
    fn default() -> Self {
        Self {
            model: WGSModel::WGS_84,
            use_afspc_mode: true,
            grav_const: get_grav_c(&WGSModel::WGS_84),
            satrec: SatRec::default(),
            use_deep_space: false,
            x2ox3: 2.0 / 3.0,
            jdstart_full: 0.0,
            jdstop_full: 0.0,
            sgp4init_out: SGP4InitOutput::default(),
        }
    }
}

impl SGP4 {
    pub fn new(model: WGSModel, use_afspc_mode: bool) -> SGP4 {
        SGP4 {
            model,
            use_afspc_mode,
            grav_const: get_grav_c(&model),
            satrec: SatRec::default(),
            use_deep_space: false,
            x2ox3: 2.0 / 3.0,
            jdstart_full: 0.0,
            jdstop_full: 0.0,
            sgp4init_out: SGP4InitOutput::default(),
        }
    }

    fn preprocess_tle(&self, tle_line1: String, tle_line2: String) -> (String, String) {
        (tle_line1, tle_line2)
    }

    /// TODO: update docstring for rust
    /// Parse TLE lines and populate SGP4 variables.
    ///
    /// This function converts the two line element (TLE) set character string data to
    /// variables and initializes the sgp4 variables. several intermediate variables
    /// and quantities are determined. The Verification mode permits quick checks of any
    /// changes to the underlying technical theory and works using a
    /// modified tle file in which the start, stop, and delta time values are
    /// included at the end of the second line of data. The Catalog mode simply
    /// propagates from -1440 to 1440 min from epoch and is useful when performing
    /// entire catalog runs.
    ///
    /// If using the FromJD mode, the start and stop Julian dates must be set before
    /// calling this function (see `set_jd_from_from_ymdhms` or `set_jd_from_yr_doy`).
    ///
    /// Args:
    ///     tle_line1 (str): First line of the TLE set
    ///     tle_line2 (str): Second line of the TLE set
    ///     typerun (TypeRun): Mode of execution (default = TypeRun.Catalog)
    ///     start (float, optional): Start time in minutes from epoch (default = None)
    ///     stop (float, optional): Stop time in minutes from epoch (default = None)
    ///     step (float, optional): Time step in minutes (default = None)
    ///
    /// Returns:
    ///     tuple (r_init, v_init, startmfe, stopmfe, deltamin)
    ///         startmfe (float): Start time in minutes from epoch
    ///         stopmfe (float): Stop time in minutes from epoch
    ///         deltamin (float): Time step in minutes
    ///         r_init (np.ndarray): Initial position vector in TEME frame in km
    ///         v_init (np.ndarray): Initial velocity vector in TEME frame in km/s
    ///
    pub fn twoline2rv(mut self, tle_line1: String, tle_line2: String) -> ([f64; 3], [f64; 3]) {
        let xpdotp = DAY2MIN / TWOPI;

        // pre-process the TLE lines
        let (tle_line1, tle_line2) = self.preprocess_tle(tle_line1, tle_line2);

        // parse the first line
        self.satrec.satnum = tle_line1[2..7].parse().unwrap();
        self.satrec.classification = Classification::from(tle_line1.chars().nth(7).unwrap());
        let binding = tle_line1[9..17].to_string();
        self.satrec.intldesg = binding.trim().to_string();
        self.satrec.epochyr = tle_line1[18..20].parse().unwrap();
        self.satrec.epochdays = tle_line1[20..32].parse().unwrap();
        self.satrec.ndot = tle_line1[33..43].parse().unwrap();
        self.satrec.nddot = f64::powi(
            tle_line1[44..50].parse().unwrap() * 10.0,
            tle_line1[50..52].parse().unwrap(),
        );
        self.satrec.bstar = f64::powi(
            tle_line1[53..59].parse().unwrap() * 10.0,
            tle_line1[59..61].parse().unwrap(),
        );
        self.satrec.elnum = tle_line1[64..68].parse().unwrap();

        // parse the second line
        self.satrec.inclo = radians(tle_line2[8..16].parse().unwrap());
        self.satrec.nodeo = radians(tle_line2[17..25].parse().unwrap());
        self.satrec.ecco = (String::from("0.") + &tle_line2[26..33].trim())
            .parse()
            .unwrap();
        self.satrec.argpo = radians(tle_line2[34..42].parse().unwrap());
        self.satrec.mo = radians(tle_line2[43..51].parse().unwrap());
        self.satrec.no_kozai = tle_line2[53..63].parse().unwrap() / xpdotp;
        self.satrec.revnum = tle_line2[63..68].parse().unwrap();

        // convert epoch year to full year
        let year = self.satrec.epochyr + if self.satrec.epochyr < 57 { 2000 } else { 1900 };

        // adjust ndot and nddot units
        self.satrec.ndot /= xpdotp * DAY2MIN;
        self.satrec.nddot /= xpdotp * DAY2MIN.powi(2);

        // compute julian date of the epoch
        let mdhms = days2mdh(year, self.satrec.epochdays);
        (self.satrec.jdsatepoch, self.satrec.jdsatepochf) = jday(year, mdhms);

        // initialize SGP4
        epoch = self.satrec.jdsatepoch + self.satrec.jdsatepochf - JD_EPOCH_1950;
        let (r_init, v_init) = self.sgp4init(epoch);

        (r_init, v_init)
    }
}
