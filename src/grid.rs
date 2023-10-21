use crate::space_domain::SpaceDomain;
use crate::cell::color;
use crate::cell::Cell;

use std::future::Future;
use std::time::{Duration, Instant};

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Text, Program};
use iced::{
    Color, Element, Length, Point, Rectangle, Renderer, Size, Theme, Vector, mouse
};
use iced::widget::canvas::event::{self, Event};


#[derive(Debug, Clone)]
pub enum Message {
    Ticked {
        tick_duration: Duration,
    },
}

#[derive(Debug, Clone)]
pub enum TickError {
    JoinFailed,
}

#[derive(Default)]
pub struct Grid {
    space_domain: SpaceDomain,
    next_cache: Cache,
    grid_cache: Cache,
    show_lines: bool,
    last_tick_duration: Duration,
    last_queued_ticks: usize,
}


impl Grid {
    pub fn are_lines_visible(&self) -> bool {
        self.show_lines
    }

    pub fn toggle_lines(&mut self, enabled: bool) {
        self.show_lines = enabled;
    }

    pub fn update(&mut self, message: Message) {
        match message {
            Message::Ticked {
                tick_duration,
            } => {
                self.next_cache.clear();

                self.last_tick_duration = tick_duration;
            }
        }
    }

    pub fn view(&self) -> Element<Message> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    pub fn tick(&mut self, amount: usize) {
        self.space_domain.tick(amount);
        self.next_cache.clear();
    }
}


impl Program<Message> for Grid {
    type State = ();


    fn draw(&self, _state: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry>{
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            self.draw_cells(frame);
        });

        if !self.show_lines {
            vec![cells]
        } else {
            let grid = self.grid_cache.draw(renderer, bounds.size(), |frame| {
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
            frame.scale(self.space_domain.get_cell_size()*2.0);

            for (x, row) in self.space_domain.get_space().iter().enumerate() {
                for (y, cell) in row.iter().enumerate() {
                    let reversed_y = row.len() - 1 - y;
                    frame.fill_rectangle(
                        Point::new(x as f32, reversed_y as f32),
                        Size::UNIT,
                        color(cell),
                    );
                    let cell_position = Point::new(x as f32, reversed_y as f32); 
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
        frame.scale(self.space_domain.get_cell_size() * 2.0);

        let total_rows = self.space_domain.get_size()[1];
        let total_columns = self.space_domain.get_size()[0];
        let width = 2.0 / self.space_domain.get_cell_size() as f32;
        let color = Color::from_rgb8(70, 74, 83);

        frame.translate(Vector::new(-width / 2.0, -width / 2.0));

        for row in 1..total_rows {
            frame.fill_rectangle(
                Point::new(0.0, row as f32),
                Size::new(total_columns as f32, width),
                color,
            );
        }

        for column in 1..total_columns {
            frame.fill_rectangle(
                Point::new(column as f32, 0.0),
                Size::new(width, total_rows as f32),
                color,
            );
        }
    }

}