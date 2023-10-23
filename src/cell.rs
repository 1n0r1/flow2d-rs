#[derive(Debug, Clone, Default)]
pub struct Cell {
    pub cell_type: CellType,
    pub velocity: [f32; 2],
    pub pressure: f32,
    pub rhs: f32,
    pub f: f32,
    pub g: f32,
    pub boundary_condition_velocity: [f32; 2]
}

#[derive(Debug, Clone, Copy, Default)]
pub enum CellType {
    #[default]
    FluidCell,
    VoidCell,
    BoundaryConditionCell(BoundaryConditionCell)
}

#[derive(Debug, Clone, Copy)]
pub enum BoundaryConditionCell {
    NoSlipCell,
    FreeSlipCell,
    OutFlowCell,
    InflowCell
}