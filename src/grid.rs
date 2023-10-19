use super::space::Space;
use super::cell::Cell;

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Text, Program};
use iced::{
    Color, Element, Length, Point, Rectangle, Renderer, Size, Theme, Vector, mouse
};
use std::time::Duration;


#[derive(Debug, Clone)]
pub enum Message {
    Ticked {
        result: Result<Space, TickError>,
        tick_duration: Duration,
    },
}

#[derive(Debug, Clone)]
pub enum TickError {
    JoinFailed,
}

#[derive(Default)]
pub struct Grid {
    space: Space,
    next_cache: Cache,
    grid_cache: Cache,
    show_lines: bool,
    last_tick_duration: Duration,
    last_queued_ticks: usize,
}



impl Program<Message> for Grid {
    type State = ();

    fn draw(&self, _state: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry>{
        let center = Vector::new(bounds.width / 2.0, bounds.height / 2.0);
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            let background = Path::rectangle(Point::ORIGIN, frame.size());
            frame.fill(&background, Color::from_rgb8(0x40, 0x44, 0x4B));

            frame.with_save(|frame| {
                frame.scale(self.space.get_cell_size());

                for (index, cell) in self.space.get_space().iter().enumerate() {
                    let column = index % self.space.get_size()[0] as usize;
                    let row = index / self.space.get_size()[1] as usize;
                
                    let cell_x = column as f32;
                    let cell_y = row as f32;
        
                    frame.fill_rectangle(
                        Point::new(cell_x, cell_y),
                        Size::UNIT,
                        Color::WHITE,
                    );
                }
            });
        });

        if self.show_lines {
            vec![cells]
        } else {
            let grid =
                self.grid_cache.draw(renderer, bounds.size(), |frame| {
                    frame.scale(self.space.get_cell_size());

                    let total_rows = self.space.get_size()[1];
                    let total_columns = self.space.get_size()[0];
                    let width = 2.0 / self.space.get_cell_size() as f32;
                    let color = Color::from_rgb8(70, 74, 83);

                    frame
                        .translate(Vector::new(-width / 2.0, -width / 2.0));

                    for row in 1..total_rows {
                        frame.fill_rectangle(
                            Point::new(0 as f32, row as f32),
                            Size::new(total_columns as f32, width),
                            color,
                        );
                    }

                    for column in 1..total_columns {
                        frame.fill_rectangle(
                            Point::new(column as f32, 0 as f32),
                            Size::new(width, total_rows as f32),
                            color,
                        );
                    }
                });
            vec![cells, grid]
        }
    }
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
                result: Ok(space),
                tick_duration,
            } => {
                // self.space.update(space);
                self.next_cache.clear();

                self.last_tick_duration = tick_duration;
            }
            Message::Ticked {
                result: Err(error), ..
            } => {
                dbg!(error);
            }
        }
    }
    
    pub fn view(&self) -> Element<Message> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}