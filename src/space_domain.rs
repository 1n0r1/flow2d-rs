use super::cell::Cell;


#[derive(Clone, Debug)]
pub struct SpaceDomain {
    space_domain: Vec<Vec<Cell>>,
    size: [u32; 2],
    cell_size: f32
}

impl Default for SpaceDomain {
    fn default() -> Self {
        let x: usize = 20;
        let y: usize = 20;
        
        let mut space_domain = vec![vec![Cell::FluidCell {
            vertical_velocity: 0.0,
            horizontal_velocity: 0.0,
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
            size: [x as u32, y as u32],
            cell_size: 20.0,
        }
    }
}

impl SpaceDomain {
    pub fn get_size(&self) -> [u32; 2] {
        self.size
    }

    pub fn get_cell_size(&self) -> f32 {
        self.cell_size
    }

    pub fn get_space(&self) -> Vec<Vec<Cell>> {
        self.space_domain.clone()
    }

    pub fn update(&self) {
        println!("updating")
    }
}