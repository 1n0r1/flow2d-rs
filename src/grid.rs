use crate::space_time_domain::SpaceTimeDomain;
use crate::cell::Cell;
use crate::cell::CellType;

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Text, Program, Stroke, Style};
use iced::{
    Color, Element, Length, Point, Rectangle, Renderer, Size, Theme, Vector, mouse
};


#[derive(Default)]
pub struct Grid {
    space_time_domain: SpaceTimeDomain,
    next_cache: Cache,
    vector_cache: Cache,
    grid_cache: Cache,
    show_lines: bool,
}


impl Grid {
    pub fn are_lines_visible(&self) -> bool {
        self.show_lines
    }

    pub fn toggle_lines(&mut self, enabled: bool) {
        self.show_lines = enabled;
    }

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
    }
}


impl Program<()> for Grid {
    type State = ();

    fn draw(&self, _state: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry>{
        const GRID_SCALE: f32 = 2000.0;
        
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            frame.scale(GRID_SCALE);
            self.draw_cells(frame);
        });
        let vectors = self.vector_cache.draw(renderer, bounds.size(), |frame| {
            frame.scale(GRID_SCALE);
            self.draw_velocity_vector(frame);
        });

        if !self.show_lines {
            vec![cells, vectors]
        } else {
            let grid = self.grid_cache.draw(renderer, bounds.size(), |frame| {
                frame.scale(GRID_SCALE);
                self.draw_grid(frame);
            });
            vec![cells, vectors, grid,]
        }
    }
}


impl Grid {
    fn draw_cells(&self, frame: &mut Frame) {
        let background = Path::rectangle(Point::ORIGIN, frame.size());
        frame.fill(&background, Color::from_rgb8(0x40, 0x44, 0x4B));

        frame.with_save(|frame| {
            let delta_x = self.space_time_domain.get_delta_space()[0];
            let delta_y = self.space_time_domain.get_delta_space()[1];
            for (x, row) in self.space_time_domain.get_space().iter().enumerate() {
                for (y, cell) in row.iter().enumerate() {
                    let pos_x = delta_x*(x as f32);
                    let reversed_y = row.len() - 1 - y;
                    let pos_y = delta_y*(reversed_y as f32);

                    frame.fill_rectangle(
                        Point::new(pos_x, pos_y),
                        Size::new(delta_x, delta_y),
                        color(cell),
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
                let pos_x = delta_x*(x as f32);
                let reversed_y = row.len() - 1 - y;
                let pos_y = delta_y*(reversed_y as f32);
                
                let velocity_scale = 0.1;
                let (u, v) = (cell.velocity[0], -cell.velocity[1]);

                let vector_start = Point::new(pos_x, pos_y);
                let vector_end = Point::new(
                    pos_x + u * velocity_scale,
                    pos_y + v * velocity_scale,
                );

                let vector = Path::line(vector_start, vector_end);
                frame.stroke(&vector, Stroke {
                    width: 1.0,
                    ..Default::default()
                });
            }
        }
    }

    fn draw_grid(&self, frame: &mut Frame) {
        let total_rows = self.space_time_domain.get_space_size()[1];
        let total_columns = self.space_time_domain.get_space_size()[0];

        let delta_x = self.space_time_domain.get_delta_space()[0];
        let delta_y = self.space_time_domain.get_delta_space()[1];
        
        let color = Color::from_rgb8(70, 74, 83);
        for row in 1..total_rows {
            let pos_row = delta_y*(row as f32);
            frame.fill_rectangle(
                Point::new(0.0, pos_row as f32),
                Size::new(total_columns as f32, 0.02),
                color,
            );
        }

        for column in 1..total_columns {
            let pos_column = delta_x*(column as f32);
            frame.fill_rectangle(
                Point::new(pos_column as f32, 0.0),
                Size::new(0.02, total_rows as f32),
                color,
            );
        }
    }

}

pub fn color(cell: &Cell) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            let hue: f32 = (0.5 + cell.pressure)*500.0 % 360.0;
            let saturation = 1.0;
            let lightness = 0.5;
            
            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)

        },
        CellType::BoundaryConditionCellType(_) => Color::from_rgb(0.5, 0.5, 0.5),
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