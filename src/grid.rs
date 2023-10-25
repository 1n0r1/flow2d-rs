use crate::simulation::Simulation;
use crate::cell;
use crate::cell::Cell;
use crate::cell::CellType;

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Program, Stroke};
use iced::{
    Color, Element, Length, Point, Rectangle, Renderer, Size, Theme, mouse
};

#[derive(Default)]
pub struct Grid {
    space_time_domain: Simulation,
    next_cache: Cache,
    vector_cache: Cache,
    contour_cache: Cache,
}


impl Grid {
    pub fn get_time(&self) -> f32{
        self.space_time_domain.get_time()
    }

    pub fn view(&self) -> Element<()> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    pub fn tick(&mut self) {
        self.space_time_domain.iterate_one_timestep();
        self.next_cache.clear();
        self.vector_cache.clear();
        self.contour_cache.clear();
    }
}


impl Program<()> for Grid {
    type State = ();

    fn draw(&self, _state: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry>{
        const GRID_SCALE: f32 = 700.0;
        
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            frame.scale(GRID_SCALE);
            self.draw_cells(frame);
        });

        let vectors = self.vector_cache.draw(renderer, bounds.size(), |frame| {
            frame.scale(GRID_SCALE);
            self.draw_velocity_vector(frame);
        });
        
        vec![cells, vectors]
    }
}


impl Grid {

    fn draw_cells(&self, frame: &mut Frame) {
        let background = Path::rectangle(Point::ORIGIN, frame.size());
        frame.fill(&background, Color::from_rgb8(0x40, 0x44, 0x4B));

        frame.with_save(|frame| {
            let delta_x = self.space_time_domain.get_delta_space()[0];
            let delta_y = self.space_time_domain.get_delta_space()[1];
            let pressure_range = self.space_time_domain.get_pressure_range();
            let speed_range = self.space_time_domain.get_speed_range();
            let psi_range = self.space_time_domain.get_psi_range();

            for (x, row) in self.space_time_domain.get_space().iter().enumerate() {
                for (y, cell) in row.iter().enumerate() {
                    let pos_x = delta_x*(x as f32);
                    let reversed_y = row.len() - 1 - y;
                    let pos_y = delta_y*(reversed_y as f32);

                    frame.fill_rectangle(
                        Point::new(pos_x, pos_y),
                        Size::new(delta_x, delta_y),
                        // color_presure(cell, pressure_range),
                        // color_speed(cell, speed_range),
                        color_psi(cell, psi_range),
                    );
                }
            }
        });
    }

    fn draw_velocity_vector(&self, frame: &mut Frame) {
        let delta_x = self.space_time_domain.get_delta_space()[0];
        let delta_y = self.space_time_domain.get_delta_space()[1];
        for (x, row) in self.space_time_domain.get_space().iter().enumerate() {
            for (y, cell) in row.iter().enumerate() {
                // if x % 2 != 0 || y % 2 != 0 {
                //     continue;
                // }
                if let cell::CellType::FluidCell = cell.cell_type {
                    let pos_x = delta_x*(x as f32);
                    let reversed_y = row.len() - 1 - y;
                    let pos_y = delta_y*(reversed_y as f32);
                    
                    let velocity_scale = 0.1;
                    let velocity = self.space_time_domain.get_centered_velocity(x, y);
                    
                    let vector_start = Point::new(pos_x + delta_x/2.0, pos_y + delta_y/2.0);
                    let vector_end = Point::new(
                        vector_start.x + velocity[0] * velocity_scale,
                        vector_start.y - velocity[1] * velocity_scale,
                    );
    
                    let vector = Path::line(vector_start, vector_end);
                    frame.stroke(&vector, Stroke {
                        width: 1.0,
                        ..Default::default()
                    });
                }
            }
        }
    }
}


pub fn color_presure(cell: &Cell, pressure_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 = 240.0 - (cell.pressure - pressure_range[0])*240.0/(pressure_range[1] - pressure_range[0]);
            let saturation = 1.0;
            let lightness = 0.5;
            
            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)

        },
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}


pub fn color_speed(cell: &Cell, speed_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            let speed = (cell.velocity[0].powi(2) + cell.velocity[1].powi(2)).sqrt();

            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 = 240.0 - (speed - speed_range[0])*240.0/(speed_range[1] - speed_range[0]);
            let saturation = 1.0;
            let lightness = 0.5;
            
            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)

        },
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}


pub fn color_psi(cell: &Cell, psi_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 = (cell.psi - psi_range[0])*1800.0/(psi_range[1] - psi_range[0])%360.0;
            let saturation = 1.0;
            let lightness = 0.5;
            
            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)

        },
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}

fn hsl_to_rgb(hue: f32, saturation: f32, lightness: f32) -> (f32, f32, f32) {
    let c = (1.0 - (2.0 * lightness - 1.0).abs()) * saturation;
    let x = c * (1.0 - ((hue / 60.0) % 2.0 - 1.0).abs());
    let m = lightness - c / 2.0;
    
    let (r, g, b) = if hue < 60.0 {
        (c, x, 0.0)
    } else if hue < 120.0 {
        (x, c, 0.0)
    } else if hue < 180.0 {
        (0.0, c, x)
    } else if hue < 240.0 {
        (0.0, x, c)
    } else if hue < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    
    (r + m, g + m, b + m)
}