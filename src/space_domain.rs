use crate::cell::BoundaryConditionCell;
use crate::cell::Cell;
use crate::cell::CellType;

pub struct SpaceDomain {
    space_domain: Vec<Vec<Cell>>,
    space_size: [usize; 2],
    delta_space: [f32; 2], // meters

    // upwind discretization parameter for evaluating spatial derivative
    gamma: f32, // 0 <= gamma <= 1

    // For coloring
    pressure_range: [f32; 2],
    speed_range: [f32; 2],
    psi_range: [f32; 2],
}

impl SpaceDomain {
    pub fn new(space_domain: Vec<Vec<Cell>>, delta_space: [f32; 2], gamma: f32) -> Self {
        let space_size = [space_domain.len(), space_domain[0].len()];
        Self {
            space_domain,
            space_size,
            delta_space,
            gamma,
            pressure_range: [0.0, 0.0],
            speed_range: [0.0, 0.0],
            psi_range: [0.0, 0.0],
        }
    }
}

// Get functions
impl SpaceDomain {
    pub fn get_delta_space(&self) -> [f32; 2] {
        self.delta_space
    }

    pub fn get_space_size(&self) -> [usize; 2] {
        self.space_size
    }

    pub fn get_space(&self) -> &Vec<Vec<Cell>> {
        &self.space_domain
    }

    pub fn get_pressure_range(&self) -> [f32; 2] {
        self.pressure_range
    }

    pub fn get_speed_range(&self) -> [f32; 2] {
        self.speed_range
    }

    pub fn get_psi_range(&self) -> [f32; 2] {
        self.psi_range
    }

    pub fn get_cell(&self, x: usize, y: usize) -> &Cell {
        &self.space_domain[x][y]
    }

    pub fn get_centered_velocity(&self, x: usize, y: usize) -> [f32; 2] {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => [
                (self.space_domain[x][y].velocity[0] + self.space_domain[x - 1][y].velocity[0])
                    / 2.0,
                (self.space_domain[x][y].velocity[1] + self.space_domain[x][y - 1].velocity[1])
                    / 2.0,
            ],
            _ => [0.0, 0.0],
        }
    }
}

// Update functions
impl SpaceDomain {
    pub fn get_cell_mut(&mut self, x: usize, y: usize) -> &mut Cell {
        &mut self.space_domain[x][y]
    }

    pub fn update_psi(&mut self) {
        self.psi_range = [0.0, 0.0];
        for x in 0..self.space_size[0] {
            self.space_domain[x][0].psi = 0.0;
            for y in 1..self.space_size[1] {
                match self.space_domain[x][y].cell_type {
                    CellType::FluidCell => {
                        self.space_domain[x][y].psi = self.space_domain[x][y - 1].psi
                            + self.space_domain[x][y].velocity[0] * self.delta_space[1];
                        if self.space_domain[x][y].psi < self.psi_range[0] {
                            self.psi_range[0] = self.space_domain[x][y].psi;
                        }
                        if self.space_domain[x][y].psi > self.psi_range[1] {
                            self.psi_range[1] = self.space_domain[x][y].psi;
                        }
                    }
                    _ => {
                        self.space_domain[x][y].psi = self.space_domain[x][y - 1].psi;
                    }
                }
            }
        }
    }

    pub fn update_pressure_and_speed_range(&mut self) {
        let (min_pressure, max_pressure, min_speed, max_speed) = self
            .space_domain
            .iter()
            .flatten()
            .filter(|cell| matches!(cell.cell_type, CellType::FluidCell))
            .map(|cell| {
                let pressure = cell.pressure;
                let speed = (cell.velocity[0].powi(2) + cell.velocity[1].powi(2)).sqrt();
                (pressure, speed)
            })
            .fold(
                (
                    f32::INFINITY,
                    f32::NEG_INFINITY,
                    f32::INFINITY,
                    f32::NEG_INFINITY,
                ),
                |(min_p, max_p, min_s, max_s), (pressure, speed)| {
                    (
                        min_p.min(pressure),
                        max_p.max(pressure),
                        min_s.min(speed),
                        max_s.max(speed),
                    )
                },
            );

        self.pressure_range = [min_pressure, max_pressure];
        self.speed_range = [min_speed, max_speed];
    }

    // Set u, v, F, G, p boundary conditions
    pub fn update_boundary_conditions(&mut self) {
        let x_size = self.space_size[0];
        let y_size = self.space_size[1];

        for x in 0..x_size {
            for y in 0..y_size {
                if let CellType::BoundaryConditionCell(bc_cell_type) =
                    self.space_domain[x][y].cell_type
                {
                    let left_cell_type: Option<CellType> =
                        (x > 0).then(|| self.space_domain[x - 1][y].cell_type);
                    let right_cell_type: Option<CellType> =
                        (x + 1 < self.space_size[0]).then(|| self.space_domain[x + 1][y].cell_type);
                    let bottom_cell_type: Option<CellType> =
                        (y > 0).then(|| self.space_domain[x][y - 1].cell_type);
                    let top_cell_type: Option<CellType> =
                        (y + 1 < self.space_size[1]).then(|| self.space_domain[x][y + 1].cell_type);

                    match bc_cell_type {
                        BoundaryConditionCell::NoSlipCell {
                            boundary_condition_velocity,
                        } => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] = 2.0
                                            * boundary_condition_velocity[1]
                                            - self.space_domain[x - 1][y].velocity[1];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x - 1][y].f =
                                            self.space_domain[x - 1][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            0.0,
                                            2.0 * boundary_condition_velocity[1]
                                                - self.space_domain[x + 1][y].velocity[1],
                                        ];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].f =
                                            self.space_domain[x][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] = 2.0
                                            * boundary_condition_velocity[0]
                                            - self.space_domain[x][y - 1].velocity[0];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y - 1].g =
                                            self.space_domain[x][y - 1].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            2.0 * boundary_condition_velocity[0]
                                                - self.space_domain[x][y + 1].velocity[0],
                                            0.0,
                                        ];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].g =
                                            self.space_domain[x][y].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure
                                    / (neighboring_fluid_count as f32)
                            }
                        }

                        BoundaryConditionCell::FreeSlipCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] =
                                            self.space_domain[x - 1][y].velocity[1];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x - 1][y].f =
                                            self.space_domain[x - 1][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity =
                                            [0.0, self.space_domain[x + 1][y].velocity[1]];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].f =
                                            self.space_domain[x][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] =
                                            self.space_domain[x][y - 1].velocity[0];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y - 1].g =
                                            self.space_domain[x][y - 1].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity =
                                            [self.space_domain[x][y + 1].velocity[0], 0.0];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].g =
                                            self.space_domain[x][y].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure
                                    / (neighboring_fluid_count as f32)
                            }
                        }

                        BoundaryConditionCell::OutFlowCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] =
                                            self.space_domain[x - 2][y].velocity[0];
                                        self.space_domain[x][y].velocity[1] =
                                            self.space_domain[x - 1][y].velocity[1];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x - 1][y].f =
                                            self.space_domain[x - 1][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            self.space_domain[x + 1][y].velocity[0],
                                            self.space_domain[x + 1][y].velocity[1],
                                        ];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].f =
                                            self.space_domain[x][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity[0] =
                                            self.space_domain[x][y - 1].velocity[0];
                                        self.space_domain[x][y - 1].velocity[1] =
                                            self.space_domain[x][y - 2].velocity[1];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y - 1].g =
                                            self.space_domain[x][y - 1].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            self.space_domain[x][y + 1].velocity[0],
                                            self.space_domain[x][y + 1].velocity[1],
                                        ];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].g =
                                            self.space_domain[x][y].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure
                                    / (neighboring_fluid_count as f32)
                            }
                        }
                        BoundaryConditionCell::InflowCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] =
                                            self.space_domain[x][y].velocity[0];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x - 1][y].f =
                                            self.space_domain[x - 1][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].f =
                                            self.space_domain[x][y].velocity[0];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] =
                                            self.space_domain[x][y].velocity[1];

                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y - 1].g =
                                            self.space_domain[x][y - 1].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].pressure +=
                                            self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].g =
                                            self.space_domain[x][y].velocity[1];
                                    }
                                    CellType::BoundaryConditionCell(_) => {}
                                    CellType::VoidCell => {}
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure
                                    / (neighboring_fluid_count as f32)
                            }
                        }
                    }
                }
            }
        }
    }
}

// Spatial derivatives
impl SpaceDomain {
    pub fn d2udx2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let ui = self.space_domain[x][y].velocity[0];
                let uip1 = self.space_domain[x + 1][y].velocity[0];
                let uim1 = self.space_domain[x - 1][y].velocity[0];
                (uip1 - 2.0 * ui + uim1) / (self.delta_space[0].powi(2))
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn d2udy2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uj = self.space_domain[x][y].velocity[0];
                let ujp1 = self.space_domain[x][y + 1].velocity[0];
                let ujm1 = self.space_domain[x][y - 1].velocity[0];
                (ujp1 - 2.0 * uj + ujm1) / (self.delta_space[1].powi(2))
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn d2vdx2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vi = self.space_domain[x][y].velocity[1];
                let vip1 = self.space_domain[x + 1][y].velocity[1];
                let vim1 = self.space_domain[x - 1][y].velocity[1];

                (vip1 - 2.0 * vi + vim1) / (self.delta_space[0].powi(2))
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn d2vdy2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vj = self.space_domain[x][y].velocity[1];
                let vjp1 = self.space_domain[x][y + 1].velocity[1];
                let vjm1 = self.space_domain[x][y - 1].velocity[1];

                (vjp1 - 2.0 * vj + vjm1) / (self.delta_space[1].powi(2))
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn du2dx(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let ui = self.space_domain[x][y].velocity[0];
                let uip1 = self.space_domain[x + 1][y].velocity[0];
                let uim1 = self.space_domain[x - 1][y].velocity[0];

                ((ui + uip1).powi(2) - (uim1 + ui).powi(2)) / 4.0 / self.delta_space[0]
                    + self.gamma
                        * ((ui + uip1).abs() * (ui - uip1) - (uim1 + ui).abs() * (uim1 - ui))
                        / 4.0
                        / self.delta_space[0]
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn dv2dy(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vj = self.space_domain[x][y].velocity[1];
                let vjp1 = self.space_domain[x][y + 1].velocity[1];
                let vjm1 = self.space_domain[x][y - 1].velocity[1];

                ((vj + vjp1).powi(2) - (vjm1 + vj).powi(2)) / 4.0 / self.delta_space[1]
                    + self.gamma
                        * ((vj + vjp1).abs() * (vj - vjp1) - (vjm1 + vj).abs() * (vjm1 - vj))
                        / 4.0
                        / self.delta_space[1]
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn duvdx(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uij = self.space_domain[x][y].velocity[0];
                let vij = self.space_domain[x][y].velocity[1];

                let vip1 = self.space_domain[x + 1][y].velocity[1];
                let vim1 = self.space_domain[x - 1][y].velocity[1];

                let uim1 = self.space_domain[x - 1][y].velocity[0];

                let ujp1 = self.space_domain[x][y + 1].velocity[0];

                let uim1jp1 = self.space_domain[x - 1][y + 1].velocity[0];

                ((uij + ujp1) * (vij + vip1) - (uim1 + uim1jp1) * (vim1 + vij))
                    / 4.0
                    / self.delta_space[0]
                    + self.gamma
                        * ((uij + ujp1).abs() * (vij - vip1)
                            - (uim1 + uim1jp1).abs() * (vim1 - vij))
                        / 4.0
                        / self.delta_space[0]
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }

    pub fn duvdy(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uij = self.space_domain[x][y].velocity[0];
                let vij = self.space_domain[x][y].velocity[1];

                let ujp1 = self.space_domain[x][y + 1].velocity[0];
                let ujm1 = self.space_domain[x][y - 1].velocity[0];

                let vjm1 = self.space_domain[x][y - 1].velocity[1];

                let vip1 = self.space_domain[x + 1][y].velocity[1];

                let vip1jm1 = self.space_domain[x + 1][y - 1].velocity[1];

                ((vij + vip1) * (uij + ujp1) - (vjm1 + vip1jm1) * (ujm1 + uij))
                    / 4.0
                    / self.delta_space[1]
                    + self.gamma
                        * ((vij + vip1).abs() * (uij - ujp1)
                            - (vjm1 + vip1jm1).abs() * (ujm1 - uij))
                        / 4.0
                        / self.delta_space[1]
            }
            _ => panic!("derivative on non fluid cell"),
        }
    }
}
