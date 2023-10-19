use super::cell::Cell;



#[derive(Clone, Debug)]
pub struct Space {
    space: Vec<Cell>,
    size: [u32; 2],
    cell_size: f32
}

impl Default for Space {
    fn default() -> Self {
        let x = 10; // Default x value
        let y = 10; // Default y value
        let total_cells = (x * y) as usize;
        let space = vec![Cell::FreeCell; total_cells];
        Self {
            space,
            size: [x, y],
            cell_size: 20.0,
        }
    }
}

impl Space {

    pub fn get_size(&self) -> [u32; 2] {
        self.size
    }

    pub fn get_cell_size(&self) -> f32 {
        self.cell_size
    }

    pub fn get_space(&self) -> Vec<Cell> {
        self.space.clone()
    }
}