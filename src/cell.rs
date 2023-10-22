use iced::Color;


#[derive(Debug, Clone)]
pub struct Cell {
    pub cell_type: CellType,
    pub velocity: [f32; 2],
    pub pressure: f32,
    pub rhs: f32,
    pub f: f32,
    pub g: f32
}

#[derive(Debug, Clone, Copy)]
pub enum CellType {
    VoidCell,
    FluidCell,
    BoundaryConditionCellType(BoundaryConditionCellType)
}

#[derive(Debug, Clone, Copy)]
pub enum BoundaryConditionCellType {
    NoSlipCell,
    FreeSlipCell,
    OutFlowCell,
    InflowCell
}

impl Cell {
    pub fn tick(&mut self, _: usize) {
        match self {
            Cell { pressure, .. } => {
                let to_add = 1.00 as f32;
                *pressure = *pressure + to_add;
            }
            _ => {}
        }
    }
}