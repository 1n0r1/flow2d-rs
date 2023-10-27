#[derive(Default)]
pub struct Cell {
    pub cell_type: CellType,
    pub velocity: [f32; 2],
    pub pressure: f32,
    pub rhs: f32,
    pub f: f32,
    pub g: f32,
    pub psi: f32,
}

#[derive(Default, Clone, Copy)]
pub enum CellType {
    #[default]
    FluidCell,
    VoidCell,
    BoundaryConditionCell(BoundaryConditionCell),
}

#[derive(Clone, Copy)]
pub enum BoundaryConditionCell {
    NoSlipCell {
        boundary_condition_velocity: [f32; 2],
    },
    FreeSlipCell,
    OutFlowCell,
    InflowCell,
}
