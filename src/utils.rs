use std::f64::consts::PI;

use crate::sgp4::{Classification, PropagationError};

pub fn radians(deg: f64) -> f64 {
    deg * PI / 180.0
}

pub fn degrees(rad: f64) -> f64 {
    rad * 180.0 / PI
}

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone)]
pub enum WGSModel {
    WGS_72_LOW_PRECISION,
    WGS_72,
    WGS_84,
}

pub struct GravitationalConstants {
    pub tumin: f64,
    pub mu: f64,
    pub radiusearthkm: f64,
    pub xke: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    pub j3oj2: f64,
}

pub fn get_grav_c(model: &WGSModel) -> GravitationalConstants {
    let (mu, rad, xke, j2, j3, j4) = match model {
        WGSModel::WGS_72 => {
            let mu = 398600.8_f64;
            let rad = 6378.135_f64;
            let xke = 60.0 / (rad.powi(3) / mu).sqrt();
            (mu, rad, xke, 0.001082616, -0.00000253881, -0.00000165597)
        }
        WGSModel::WGS_84 => {
            let mu = 398600.5_f64;
            let rad = 6378.137_f64;
            let xke = 60.0 / (rad.powi(3) / mu).sqrt();
            (
                mu,
                rad,
                xke,
                0.00108262998905,
                -0.00000253215306,
                -0.00000161098761,
            )
        }
        WGSModel::WGS_72_LOW_PRECISION => (
            398600.79964,
            6378.135,
            0.0743669161,
            0.001082616,
            -0.00000253881,
            -0.00000165597,
        ),
    };
    GravitationalConstants {
        tumin: 1.0 / xke,
        mu,
        radiusearthkm: rad,
        xke,
        j2,
        j3,
        j4,
        j3oj2: j3 / j2,
    }
}

#[derive(Default)]
pub struct SatRec {
    pub a: f64,
    pub alta: f64,
    pub altp: f64,
    pub mo: f64,
    pub mdot: f64,
    pub argpo: f64,
    pub argpdot: f64,
    pub nodeo: f64,
    pub nodedot: f64,
    pub nodecf: f64,
    pub cc1: f64,
    pub cc4: f64,
    pub cc5: f64,
    pub omgcof: f64,
    pub xmcod: f64,
    pub eta: f64,
    pub sinmao: f64,
    pub delmo: f64,
    pub d2: f64,
    pub d3: f64,
    pub d4: f64,
    pub t2cof: f64,
    pub t3cof: f64,
    pub t4cof: f64,
    pub t5cof: f64,
    pub no: f64,
    pub ecco: f64,
    pub inclo: f64,
    pub isimp: bool,
    pub bstar: f64,
    pub xfact: f64,
    pub xlamo: f64,
    pub atime: f64,
    pub error: Option<PropagationError>,
    pub t: f64,
    pub aycof: f64,
    pub xlcof: f64,
    pub x1mth2: f64,
    pub x7thm1: f64,
    pub satnum: u64,
    pub intldesg: String,
    pub epochyr: u64,
    pub epochdays: f64,
    pub ndot: f64,
    pub nddot: f64,
    pub elnum: u64,
    pub revnum: u64,
    pub no_kozai: f64,
    pub jdsatepoch: f64,
    pub jdsatepochf: f64,
    pub init: bool,
    pub classification: Option<Classification>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_radians() {
        assert_eq!(radians(0.0), 0.0);
        assert_eq!(radians(90.0), PI / 2.0);
        assert_eq!(radians(45.0), PI / 4.0);
        assert_eq!(radians(270.0), 3.0 * PI / 2.0);
    }

    #[test]
    fn test_degrees() {
        assert_eq!(degrees(0.0), 0.0);
        assert_eq!(degrees(PI / 2.0), 90.0);
        assert_eq!(degrees(3.0 * PI / 2.0), 270.0);
    }
}
