use crate::cell::Cell;
use crate::cell::CellType;
use crate::cell::BoundaryConditionCell;


pub struct SpaceTimeDomainPreset {
    pub space_domain: Vec<Vec<Cell>>,
    pub space_size: [usize; 2],
    pub delta_space: [f32; 2], // meters
    pub delta_time: f32, // seconds,
    pub acceleration: [f32; 2], // meters/seconds^2
    pub reynolds: f32,
}

pub fn lid_driven_cavity() -> SpaceTimeDomainPreset {
    let x_length = 1.0;
    let y_length = 1.0;
    let x: usize = 128;
    let y: usize = 128;
    
    let mut space_domain = vec![vec![Cell::default(); y]; x];
    for xi in 0..x {
        for yi in 0..y {
            if xi == 0 || xi == x - 1 || yi == 0 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    ..Default::default()
                };
            }
            if yi == y - 1 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    boundary_condition_velocity: [1.0, 0.0],
                    ..Default::default()
                }
            }
        }
    }
    for xi in [0, x - 1] {
        for yi in [0, y - 1] {
            space_domain[xi][yi] = Cell {
                cell_type: CellType::VoidCell,
                ..Default::default()
            };

        }
    }
    SpaceTimeDomainPreset {
        space_domain,
        space_size: [x, y],
        delta_space: [x_length/(x as f32), y_length/(y as f32)],
        delta_time: 0.005,
        reynolds: 1000.0,
        acceleration: [0.0, 0.0]
    }
}


pub fn backward_facing_step() -> SpaceTimeDomainPreset {
    let x_length = 15.0;
    let y_length = 1.5;
    let x: usize = 150;
    let y: usize = 75;

    let mut space_domain = vec![vec![Cell::default(); y]; x];
    for xi in 0..x {
        for yi in 0..y {
            if xi == 0 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::InflowCell),
                    velocity: [1.0, 0.0],
                    ..Default::default()
                };
                continue;
            }
            if xi == x - 1 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::OutFlowCell),
                    ..Default::default()
                };
                continue;
            }
            if yi == y - 1 || yi == 0 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    ..Default::default()
                };
                continue;
            }
        }
    }
    for xi in [0, x - 1] {
        for yi in [0, y - 1] {
            space_domain[xi][yi] = Cell {
                cell_type: CellType::VoidCell,
                ..Default::default()
            };

        }
    }

    for xi in 0..75 {
        for yi in 0..37 {
            space_domain[xi][yi] = Cell {
                cell_type: CellType::VoidCell,
                ..Default::default()
            };
        }
    }

    for xi in 0..76 {
        for yi in 0..38 {
            if xi == 75 || yi == 37{
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    ..Default::default()
                };
            }
        }
    }

    SpaceTimeDomainPreset {
        space_domain,
        space_size: [x, y],
        delta_space: [x_length/(x as f32), y_length/(y as f32)],
        delta_time: 0.005,
        reynolds: 1000.0,
        acceleration: [0.0, 0.0]
    }
}


pub fn cylinder_cross_flow() -> SpaceTimeDomainPreset {
    let x_length = 11.0;
    let y_length = 4.1;
    let x: usize = 110;
    let y: usize = 41;

    let mut space_domain = vec![vec![Cell::default(); y]; x];
    for xi in 0..x {
        for yi in 0..y {
            if xi == 0 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::InflowCell),
                    velocity: [1.5, 0.0],
                    ..Default::default()
                };
                continue;
            }
            if xi == x - 1 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::OutFlowCell),
                    ..Default::default()
                };
                continue;
            }
            if yi == y - 1 || yi == 0 {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    ..Default::default()
                };
                continue;
            }
        }
    }
    for xi in [0, x - 1] {
        for yi in [0, y - 1] {
            space_domain[xi][yi] = Cell {
                cell_type: CellType::VoidCell,
                ..Default::default()
            };

        }
    }

    let delta_space = [x_length/(x as f32), y_length/(y as f32)];

    let radius = 5.0;
    let center = [20, 20];
    
    for xi in 14..26 as i32 {
        for yi in 14..26 as i32 {
            let x_dist = xi - center[0];
            let y_dist = yi - center[1];
            let distance = ((x_dist * x_dist + y_dist * y_dist) as f32 ).sqrt();
    
            if distance <= radius {
                space_domain[xi as usize][yi as usize] = Cell {
                    cell_type: CellType::BoundaryConditionCell(BoundaryConditionCell::NoSlipCell),
                    ..Default::default()
                };
            }
        }
    }
    
    
    SpaceTimeDomainPreset {
        space_domain,
        space_size: [x, y],
        delta_space: delta_space,
        delta_time: 0.001,
        reynolds: 100.0,
        acceleration: [0.0, 0.0]
    }
}
