use iced::Color;


#[derive(Debug, Clone)]
pub struct Cell {
    pub cell_type: CellType,
    pub velocity: [f32; 2],
    pub pressure: f32,
    pub rhs: f32,
    pub f: f32,
    pub g: f32,
    pub boundary_condition_velocity: [f32; 2]
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