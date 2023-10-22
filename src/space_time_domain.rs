use crate::cell::Cell;
use crate::cell::CellType;
use crate::cell::BoundaryConditionCellType;

#[derive(Debug, Clone)]
pub struct SpaceTimeDomain {
    space_domain: Vec<Vec<Cell>>,
    space_size: [usize; 2],
    time: f32, // seconds
    delta_space: [f32; 2], // meters
    delta_time: f32, // seconds,
    reynolds: f32,
    acceleration: [f32; 2], // meters/seconds^2
    gamma: f32
}

impl Default for SpaceTimeDomain {
    fn default() -> Self {
        let x: usize = 20;
        let y: usize = 20;
        
        let mut space_domain = vec![vec![Cell {
            cell_type: CellType::FluidCell,
            velocity: [0.0, 0.0],
            pressure: 0.0,
            rhs: 0.0,
            f: 0.0,
            g: 0.0,
        }; y]; x];
        
        for xi in 0..x {
            for yi in 0..y {
                if xi == 0 || xi == x - 1 || yi == 0 || yi == y - 1 {
                    space_domain[xi][yi] = Cell {
                        cell_type: CellType::BoundaryConditionCellType(BoundaryConditionCellType::NoSlipCell),
                        velocity: [0.0, 0.0],
                        pressure: 0.0,
                        rhs: 0.0,
                        f: 0.0,
                        g: 0.0,
                    };
                }
            }
        }

        for xi in [0, x - 1] {
            for yi in [0, y - 1] {
                space_domain[xi][yi] = Cell {
                    cell_type: CellType::VoidCell,
                    velocity: [0.0, 0.0],
                    pressure: 0.0,
                    rhs: 0.0,
                    f: 0.0,
                    g: 0.0,
                };

            }
        }


        Self {
            space_domain,
            space_size: [x, y],
            time: 0.0,
            delta_space: [1.0, 1.0],
            delta_time: 0.1,
            reynolds: 100.0,
            acceleration: [0.0, -9.8],
            gamma: 0.1
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

    pub fn tick(&mut self, amount: usize) {
        self.set_boundary_conditions();

        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                
                self.space_domain[x as usize][y as usize].tick(amount);
            }
        }
        self.time += self.delta_time*(amount as f32)
    }

}


impl SpaceTimeDomain {
    fn calculate_and_update_fg(&mut self, amount: usize) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                match self.space_domain[x][y].cell_type {
                    CellType::FluidCell => {
                        
                        
                    },
                    _ => {}
                }
            }
        }
        self.time += self.delta_time*(amount as f32)
    }
    

    fn set_boundary_conditions(&mut self) {
        let x_size = self.space_size[0];
        let y_size = self.space_size[1];

        for x in 0..x_size {
            for y in 0..y_size {
                let cell_type = &self.space_domain[x][y].cell_type;

                let left_cell_type: Option<CellType> = (x > 0).then(|| self.space_domain[x - 1][y].cell_type);
                let right_cell_type: Option<CellType> = (x + 1 < self.space_size[0]).then(|| self.space_domain[x + 1][y].cell_type);
                let bottom_cell_type: Option<CellType> = (y > 0).then(|| self.space_domain[x][y - 1].cell_type);
                let top_cell_type: Option<CellType> = (y + 1 < self.space_size[1]).then(|| self.space_domain[x][y + 1].cell_type);

                if let CellType::BoundaryConditionCellType(bc_cell_type) = cell_type {
                    match bc_cell_type {
                        BoundaryConditionCellType::NoSlipCell => {
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] = -self.space_domain[x - 1][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [0.0, -self.space_domain[x + 1][y].velocity[1]];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] = -self.space_domain[x][y - 1].velocity[0];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [-self.space_domain[x][y + 1].velocity[0], 0.0]
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                        },

                        BoundaryConditionCellType::FreeSlipCell => {
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = 0.0;
                                        self.space_domain[x][y].velocity[1] = self.space_domain[x - 1][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [0.0, self.space_domain[x + 1][y].velocity[1]];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y - 1].velocity[1] = 0.0;
                                        self.space_domain[x][y].velocity[0] = self.space_domain[x][y - 1].velocity[0];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x][y + 1].velocity[0], 0.0]
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                        },

                        BoundaryConditionCellType::OutFlowCell => {
                            if let Some(left_cell_type) = left_cell_type {
                                match left_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x - 1][y].velocity[0] = self.space_domain[x - 2][y].velocity[0];
                                        self.space_domain[x][y].velocity[1] = self.space_domain[x - 1][y].velocity[1];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(right_cell_type) = right_cell_type {
                                match right_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x + 1][y].velocity[0], self.space_domain[x + 1][y].velocity[1]];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(bottom_cell_type) = bottom_cell_type {
                                match bottom_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity[0] = self.space_domain[x][y - 1].velocity[0];
                                        self.space_domain[x][y - 1].velocity[1] = self.space_domain[x][y - 2].velocity[1];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                            if let Some(top_cell_type) = top_cell_type {
                                match top_cell_type {
                                    CellType::FluidCell => {
                                        self.space_domain[x][y].velocity = [self.space_domain[x][y + 1].velocity[0], self.space_domain[x][y + 1].velocity[1]];
                                    },
                                    CellType::BoundaryConditionCellType(_) => {},
                                    CellType::VoidCell => {},
                                }
                            }
                        },
                        BoundaryConditionCellType::InflowCell => {
                            todo!()
                        },
                    }
                }
                
            }
        }
    }

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
                self.gamma*((ui + uip1).abs()*(ui-uip1) - (uim1 + ui).abs()*(uim1-ui))/4.0/self.delta_space[0]

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
                    self.gamma*((vj + vjp1).abs()*(vj-vjp1) - (vjm1 + vj).abs()*(vjm1-vj))/4.0/self.delta_space[1]

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
                    self.gamma*((uij + ujp1).abs()*(vij - vip1) - (uim1 + uim1jp1).abs()*(vim1 - vij))/4.0/self.delta_space[0]

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
                    self.gamma*((vij + vip1).abs()*(uij - ujp1) - (vjm1 + vip1jm1).abs()*(ujm1 - uij))/4.0/self.delta_space[1]

            },
            _ => panic!("derivative on non fluid cell"),
        }
    }

}