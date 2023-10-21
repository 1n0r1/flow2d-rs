use crate::space_time_domain::SpaceTimeDomain;
use crate::cell::color;
use crate::cell::Cell;

use std::future::Future;
use std::time::{Duration, Instant};

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Text, Program};
use iced::{
    Color, Element, Length, Point, Rectangle, Renderer, Size, Theme, Vector, mouse
};


#[derive(Default)]
pub struct Grid {
    space_time_domain: SpaceTimeDomain,
    next_cache: Cache,
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

    pub fn view(&self) -> Element<()> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    pub fn tick(&mut self, amount: usize) {
        self.space_time_domain.tick(amount);
        self.next_cache.clear();
    }
}


impl Program<()> for Grid {
    type State = ();


    fn draw(&self, _state: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry>{
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            frame.scale(40.0);
            self.draw_cells(frame);
        });

        if !self.show_lines {
            vec![cells]
        } else {
            let grid = self.grid_cache.draw(renderer, bounds.size(), |frame| {
                frame.scale(40.0);
                self.draw_grid(frame);
            });
            vec![cells, grid]
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
                    let cell_position = Point::new(pos_x, pos_y); 
                    let label_position = cell_position + Vector::new(0.1, 0.1);
                    let label_text = format!("({}, {})", x, y);
                    let label = Text {
                        content: label_text,
                        position: label_position,
                        color: Color::BLACK,
                        size: 10.0,
                        ..Text::default()
                    };
                    
                    frame.fill_text(label);
                }
            }
        });
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