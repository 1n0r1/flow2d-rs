use crate::cell::Cell;
use crate::cell::CellType;
use crate::cell::BoundaryConditionCell;

const GAMMA: f32 = 0.9; // 0 <= GAMMA <= 1
const OMEGA: f32 = 1.7; // 0 <= OMEGA <= 2
const ITR_MAX: usize = 100;
const POISSON_EPSILON: f32 = 0.001;

pub struct SpaceTimeDomain {
    space_domain: Vec<Vec<Cell>>,
    space_size: [usize; 2],
    time: f32, // seconds
    delta_space: [f32; 2], // meters
    delta_time: f32, // seconds,
    acceleration: [f32; 2], // meters/seconds^2
    reynolds: f32,
}
impl Default for SpaceTimeDomain {
    fn default() -> Self {
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


        Self {
            space_domain,
            space_size: [x, y],
            time: 0.0,
            delta_space: [1.0/(x as f32), 1.0/(y as f32)],
            delta_time: 0.005,
            reynolds: 1000.0,
            acceleration: [0.0, 0.0]
        }
    }
}

impl SpaceTimeDomain {
    pub fn get_space_size(&self) -> [usize; 2] {
        self.space_size
    }

    pub fn get_delta_space(&self) -> [f32; 2] {
        self.delta_space
    }

    pub fn get_space(&self) -> Vec<Vec<Cell>> {
        self.space_domain.clone()
    }

    pub fn get_time(&self) -> f32 {
        self.time
    }
    
    pub fn iterate_one_timestep(&mut self) {
        self.set_boundary_conditions();
        self.update_fg();
        self.update_rhs();
        self.solve_poisson_pressure_equation();
        self.update_velocity();

        self.time += self.delta_time
    }

}


impl SpaceTimeDomain {
    fn update_velocity(&mut self) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                match self.space_domain[x][y].cell_type {
                    CellType::FluidCell => {
                        let right_cell_type: Option<CellType> = (x + 1 < self.space_size[0]).then(|| self.space_domain[x + 1][y].cell_type);
                        let top_cell_type: Option<CellType> = (y + 1 < self.space_size[1]).then(|| self.space_domain[x][y + 1].cell_type);
        
                        if let Some(right_cell_type) = right_cell_type {
                            if let CellType::BoundaryConditionCell(_) = right_cell_type {
                            
                            } else {
                                self.space_domain[x][y].velocity[0] = self.space_domain[x][y].f - self.delta_time*(self.space_domain[x + 1][y].pressure - self.space_domain[x][y].pressure)/self.delta_space[0]
                            }
                        }

                        if let Some(top_cell_type) = top_cell_type {
                            if let CellType::BoundaryConditionCell(_) = top_cell_type {
                            
                            } else {
                                self.space_domain[x][y].velocity[1] = self.space_domain[x][y].g - self.delta_time*(self.space_domain[x][y + 1].pressure - self.space_domain[x][y].pressure)/self.delta_space[1]
                            }
                            
                        }
                    },
                    _ => {}
                }

            }
        }
    }
    
    fn solve_poisson_pressure_equation(&mut self) {
        for _ in 0..ITR_MAX {
            let mut residual_norm : f32 = 0.0;
            
            for x in 0..self.space_size[0] {
                for y in 0..self.space_size[1] {
                    match self.space_domain[x][y].cell_type {
                        CellType::FluidCell => {
                            residual_norm += (
                                (self.space_domain[x + 1][y].pressure - 2.0*self.space_domain[x][y].pressure + self.space_domain[x - 1][y].pressure 
                                )/self.delta_space[0].powi(2)
                                + (self.space_domain[x][y + 1].pressure - 2.0*self.space_domain[x][y].pressure + self.space_domain[x][y - 1].pressure 
                                )/self.delta_space[1].powi(2)
                                - self.space_domain[x][y].rhs
                            ).powi(2);
                        },
                        _ => {}
                    }

                }
            }
            if (residual_norm/(self.space_size[0] as f32)/(self.space_size[1] as f32)).sqrt() < POISSON_EPSILON {
                break;
            }

            self.update_pressures_for_boundary_cells();

            for x in 0..self.space_size[0] {
                for y in 0..self.space_size[1] {
                    match self.space_domain[x][y].cell_type {
                        CellType::FluidCell => {
                            self.space_domain[x][y].pressure = 
                                (1.0 - OMEGA)*self.space_domain[x][y].pressure + OMEGA*(
                                    (self.space_domain[x + 1][y].pressure + (self.space_domain[x - 1][y].pressure))/self.delta_space[0].powi(2)
                                    + (self.space_domain[x][y + 1].pressure + (self.space_domain[x][y - 1].pressure))/self.delta_space[1].powi(2)
                                    - self.space_domain[x][y].rhs
                                )
                                /(2.0/self.delta_space[0].powi(2) + 2.0/self.delta_space[1].powi(2))
                            ;
                        },
                        _ => {}
                    }
                }
            }


        }
    }

    fn update_pressures_for_boundary_cells(&mut self) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                let cell_type = &self.space_domain[x][y].cell_type;

                if let CellType::BoundaryConditionCell(_) = cell_type {
                    let neighboring_cells = [
                        (x.wrapping_sub(1), y),
                        (x + 1, y),
                        (x, y.wrapping_sub(1)),
                        (x, y + 1),
                    ];
                    
                    let mut neighboring_fluid_count = 0;
                    self.space_domain[x][y].pressure = 0.0;

                    for (dx, dy) in neighboring_cells.iter() {
                        if let Some(cell) = self.space_domain.get(*dx).and_then(|row| row.get(*dy)) {
                            match cell.cell_type {
                                CellType::FluidCell => {
                                    self.space_domain[x][y].pressure += cell.pressure;
                                    neighboring_fluid_count += 1;
                                }
                                _ => {}
                            }
                        }
                    }
                    if neighboring_fluid_count != 0 {
                        self.space_domain[x][y].pressure = self.space_domain[x][y].pressure/(neighboring_fluid_count as f32);
                    }
                }
            }
        }
    }

    fn update_rhs(&mut self) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                match self.space_domain[x][y].cell_type {
                    CellType::FluidCell => {
                        self.space_domain[x][y].rhs = (
                            (self.space_domain[x][y].f - self.space_domain[x - 1][y].f)/self.delta_space[0]
                            + (self.space_domain[x][y].g - self.space_domain[x][y - 1].g)/self.delta_space[1]
                        )/self.delta_time;
                    },
                    _ => { }
                }
            }
        }
    }

    fn update_fg(&mut self) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                match self.space_domain[x][y].cell_type {
                    CellType::FluidCell => {
                        let right_cell_type: Option<CellType> = (x + 1 < self.space_size[0]).then(|| self.space_domain[x + 1][y].cell_type);
                        let top_cell_type: Option<CellType> = (y + 1 < self.space_size[1]).then(|| self.space_domain[x][y + 1].cell_type);
        
                        if let Some(right_cell_type) = right_cell_type {
                            if let CellType::BoundaryConditionCell(_) = right_cell_type {
                            
                            } else {
                                self.space_domain[x][y].f = self.space_domain[x][y].velocity[0]
                                + self.delta_time*(
                                    (self.d2udx2(x, y) + self.d2udy2(x, y))/self.reynolds
                                    - self.du2dx(x, y) - self.duvdy(x, y) + self.acceleration[0]
                                );
                            }
                        }

                        if let Some(top_cell_type) = top_cell_type {
                            if let CellType::BoundaryConditionCell(_) = top_cell_type {
                            
                            } else {
                                self.space_domain[x][y].g= self.space_domain[x][y].velocity[1]
                                + self.delta_time*(
                                    (self.d2vdx2(x, y) + self.d2vdy2(x, y))/self.reynolds
                                    - self.duvdx(x, y) - self.dv2dy(x, y) + self.acceleration[1]
                                )
                            }
                            
                        }

                    },
                    _ => { }
                }
            }
        }
    }
    

    fn set_boundary_conditions(&mut self) {
        // Set u, v, F, G, p boundary conditions
        
        let x_size = self.space_size[0];
        let y_size = self.space_size[1];

        for x in 0..x_size {
            for y in 0..y_size {
                let cell_type = &self.space_domain[x][y].cell_type;

                if let CellType::BoundaryConditionCell(bc_cell_type) = cell_type {
                    let left_cell_type: Option<CellType> = (x > 0).then(|| self.space_domain[x - 1][y].cell_type);
                    let right_cell_type: Option<CellType> = (x + 1 < self.space_size[0]).then(|| self.space_domain[x + 1][y].cell_type);
                    let bottom_cell_type: Option<CellType> = (y > 0).then(|| self.space_domain[x][y - 1].cell_type);
                    let top_cell_type: Option<CellType> = (y + 1 < self.space_size[1]).then(|| self.space_domain[x][y + 1].cell_type);

                    match bc_cell_type {
                        BoundaryConditionCell::NoSlipCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] = 
                                            2.0*self.space_domain[x][y].boundary_condition_velocity[1]
                                            - self.space_domain[x - 1][y].velocity[1];

                                        self.space_domain[x][y].pressure += self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x - 1][y].f = self.space_domain[x - 1][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            0.0, 
                                            2.0*self.space_domain[x][y].boundary_condition_velocity[1] - self.space_domain[x + 1][y].velocity[1]
                                        ];

                                        self.space_domain[x][y].pressure += self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;

                                        self.space_domain[x][y].f = self.space_domain[x][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] = 
                                            2.0*self.space_domain[x][y].boundary_condition_velocity[0]
                                            - self.space_domain[x][y - 1].velocity[0];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y - 1].g = self.space_domain[x][y - 1].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [
                                            2.0*self.space_domain[x][y].boundary_condition_velocity[0] - self.space_domain[x][y + 1].velocity[0],
                                            0.0
                                        ];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y].g = self.space_domain[x][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure/(neighboring_fluid_count as f32)
                            }
                        },

                        BoundaryConditionCell::FreeSlipCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] = self.space_domain[x - 1][y].velocity[1];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x - 1][y].f = self.space_domain[x - 1][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [0.0, self.space_domain[x + 1][y].velocity[1]];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y].f = self.space_domain[x][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] = self.space_domain[x][y - 1].velocity[0];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y - 1].g = self.space_domain[x][y - 1].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x][y + 1].velocity[0], 0.0];

                                        self.space_domain[x][y].pressure += self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y].g = self.space_domain[x][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure/(neighboring_fluid_count as f32)
                            }
                        },

                        BoundaryConditionCell::OutFlowCell => {
                            self.space_domain[x][y].pressure = 0.0;
                            let mut neighboring_fluid_count = 0;
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = self.space_domain[x - 2][y].velocity[0];
                                        self.space_domain[x][y].velocity[1] = self.space_domain[x - 1][y].velocity[1];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x - 1][y].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x - 1][y].f = self.space_domain[x - 1][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x + 1][y].velocity[0], self.space_domain[x + 1][y].velocity[1]];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x + 1][y].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y].f = self.space_domain[x][y].velocity[0];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity[0] = self.space_domain[x][y - 1].velocity[0];
                                        self.space_domain[x][y - 1].velocity[1] = self.space_domain[x][y - 2].velocity[1];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x][y - 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y - 1].g = self.space_domain[x][y - 1].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x][y + 1].velocity[0], self.space_domain[x][y + 1].velocity[1]];
                                        
                                        self.space_domain[x][y].pressure += self.space_domain[x][y + 1].pressure;
                                        neighboring_fluid_count += 1;
                                        
                                        self.space_domain[x][y].g = self.space_domain[x][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCell(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if neighboring_fluid_count != 0 {
                                self.space_domain[x][y].pressure = self.space_domain[x][y].pressure/(neighboring_fluid_count as f32)
                            }
                        },
                        BoundaryConditionCell::InflowCell => {
                            todo!()
                        },
                    }
                }
                
            }
        }
    }

}

// Derivative stuffs
impl SpaceTimeDomain {
    fn d2udx2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let ui = self.space_domain[x][y].velocity[0];
                let uip1 = self.space_domain[x + 1][y].velocity[0];
                let uim1 = self.space_domain[x - 1][y].velocity[0];
                (uip1 - 2.0*ui + uim1)/(self.delta_space[0].powi(2))
            }
            _ => panic!("derivative on non fluid cell")
        }
    }
    
    fn d2udy2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uj = self.space_domain[x][y].velocity[0];
                let ujp1 = self.space_domain[x][y + 1].velocity[0];
                let ujm1 = self.space_domain[x][y - 1].velocity[0];
                (ujp1 - 2.0*uj + ujm1)/(self.delta_space[1].powi(2))
            }
            _ => panic!("derivative on non fluid cell")
        }
    }
    
    fn d2vdx2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vi = self.space_domain[x][y].velocity[1];
                let vip1 = self.space_domain[x + 1][y].velocity[1];
                let vim1 = self.space_domain[x - 1][y].velocity[1];

                (vip1 - 2.0*vi + vim1)/(self.delta_space[0].powi(2))
            },
            _ => panic!("derivative on non fluid cell"),
        }
    }
    
    fn d2vdy2(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vj = self.space_domain[x][y].velocity[1];
                let vjp1 = self.space_domain[x][y + 1].velocity[1];
                let vjm1 = self.space_domain[x][y - 1].velocity[1];
                
                (vjp1 - 2.0*vj + vjm1)/(self.delta_space[1].powi(2))
            },
            _ => panic!("derivative on non fluid cell"),
        }
    }

    fn du2dx(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let ui = self.space_domain[x][y].velocity[0];
                let uip1 = self.space_domain[x + 1][y].velocity[0];
                let uim1 = self.space_domain[x - 1][y].velocity[0];
                
                ((ui + uip1).powi(2) - (uim1 + ui).powi(2))/4.0/self.delta_space[0] +
                GAMMA*((ui + uip1).abs()*(ui-uip1) - (uim1 + ui).abs()*(uim1-ui))/4.0/self.delta_space[0]

            },
            _ => panic!("derivative on non fluid cell"),
        }
    }
    
    fn dv2dy(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let vj = self.space_domain[x][y].velocity[1];
                let vjp1 = self.space_domain[x][y + 1].velocity[1];
                let vjm1 = self.space_domain[x][y - 1].velocity[1];
                
                ((vj + vjp1).powi(2) - (vjm1 + vj).powi(2))/4.0/self.delta_space[1] +
                    GAMMA*((vj + vjp1).abs()*(vj-vjp1) - (vjm1 + vj).abs()*(vjm1-vj))/4.0/self.delta_space[1]

            },
            _ => panic!("derivative on non fluid cell"),
        }
    }
    
    fn duvdx(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uij = self.space_domain[x][y].velocity[0];
                let vij = self.space_domain[x][y].velocity[1];

                let vip1 = self.space_domain[x + 1][y].velocity[1];
                let vim1 = self.space_domain[x - 1][y].velocity[1];

                let uim1 = self.space_domain[x - 1][y].velocity[0];
                
                let ujp1 = self.space_domain[x][y + 1].velocity[0];
                
                let uim1jp1 = self.space_domain[x - 1][y + 1].velocity[0];

                ((uij + ujp1)*(vij + vip1) - (uim1 + uim1jp1)*(vim1 + vij))/4.0/self.delta_space[0] +
                    GAMMA*((uij + ujp1).abs()*(vij - vip1) - (uim1 + uim1jp1).abs()*(vim1 - vij))/4.0/self.delta_space[0]

            },
            _ => panic!("derivative on non fluid cell"),
        }
    }

    fn duvdy(&self, x: usize, y: usize) -> f32 {
        match self.space_domain[x][y].cell_type {
            CellType::FluidCell => {
                let uij = self.space_domain[x][y].velocity[0];
                let vij = self.space_domain[x][y].velocity[1];
                
                let ujp1 = self.space_domain[x][y + 1].velocity[0];
                let ujm1 = self.space_domain[x][y - 1].velocity[0];

                let vjm1 = self.space_domain[x][y - 1].velocity[1];

                let vip1 = self.space_domain[x + 1][y].velocity[1];
                
                let vip1jm1 = self.space_domain[x + 1][y - 1].velocity[1];

                ((vij + vip1)*(uij + ujp1) - (vjm1 + vip1jm1)*(ujm1 + uij))/4.0/self.delta_space[1] +
                    GAMMA*((vij + vip1).abs()*(uij - ujp1) - (vjm1 + vip1jm1).abs()*(ujm1 - uij))/4.0/self.delta_space[1]

            },
            _ => panic!("derivative on non fluid cell"),
        }
    }
}