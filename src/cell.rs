

#[derive(Clone, Debug)]
pub enum Cell {
    FreeCell,
    WaterCell,
    BoundaryConditionCell,
}

pub struct WaterCell {
    vertical_velocity: f64,
    horizontal_velocity: f64,
    pressure: f64,
    temperature: f64
}