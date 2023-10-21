use crate::cell::Cell;


#[derive(Debug, Clone)]
pub struct SpaceTimeDomain {
    space_domain: Vec<Vec<Cell>>,
    space_size: [u32; 2],
    delta_space: [f32; 2], // meters
    delta_time: f32, // seconds
}

impl Default for SpaceTimeDomain {
    fn default() -> Self {
        let x: usize = 20;
        let y: usize = 20;
        
        let mut space_domain = vec![vec![Cell::FluidCell {
            velocity: [0.0, 0.0],
            pressure: 0.0,
            temperature: 0.0,
            rhs: 0.0,
            f: 0.0,
            g: 0.0,
        }; y]; x];
        
        for xi in 0..x {
            for yi in 0..y {
                if xi == 0 || xi == x - 1 || yi == 0 || yi == y - 1 {
                    space_domain[xi][yi] = Cell::BoundaryConditionCell;
                }
            }
        }

        Self {
            space_domain,
            space_size: [x as u32, y as u32],
            delta_space: [1.0, 1.0],
            delta_time: 0.1
        }
    }
}

impl SpaceTimeDomain {
    pub fn get_space_size(&self) -> [u32; 2] {
        self.space_size
    }

    pub fn get_delta_space(&self) -> [f32; 2] {
        self.delta_space
    }

    pub fn get_space(&self) -> Vec<Vec<Cell>> {
        self.space_domain.clone()
    }

    pub fn tick(&mut self, amount: usize) {
        for x in 0..self.space_size[0] {
            for y in 0..self.space_size[1] {
                self.space_domain[x as usize][y as usize].tick(amount);
            }
        }
    }

}