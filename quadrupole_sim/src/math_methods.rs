use nalgebra::{Vector2, vector};

pub fn x_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],      // x'
        -k * state[0], // x''
    ]
}

pub fn y_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],     // y'
        k * state[0], // y''
    ]
}

pub fn rk4_step<F>(state: Vector2<f64>, z: f64, dz: f64, f: F) -> Vector2<f64>
where
    F: Fn(f64, Vector2<f64>) -> Vector2<f64>,
{
    let k1 = f(z, state);
    let k2 = f(z + 0.5 * dz, state + 0.5 * dz * k1);
    let k3 = f(z + 0.5 * dz, state + 0.5 * dz * k2);
    let k4 = f(z + dz, state + dz * k3);

    state + (dz / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}
