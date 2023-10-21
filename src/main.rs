mod grid;
mod space_domain;
mod cell;

use grid::Grid;


use iced::executor;
use iced::theme::{self, Theme};
use iced::time;
use iced::widget::{
    button, checkbox, column, container, row, slider, text,
};
use iced::window;
use iced::{
    Alignment, Application, Command, Element, Length, Settings, Subscription,
};
use std::time::{Duration, Instant};

pub fn main() -> iced::Result {
    EulerFluidSimulation::run(Settings {
        antialiasing: true,
        window: window::Settings {
            position: window::Position::Centered,
            ..window::Settings::default()
        },
        ..Settings::default()
    })
}

#[derive(Default)]
struct EulerFluidSimulation {
    grid: Grid,
    is_playing: bool,
    queued_ticks: usize,
    speed: usize,
    next_speed: Option<usize>,
    version: usize,
}

#[derive(Debug, Clone)]
enum Message {
    Grid(grid::Message, usize),
    Tick(Instant),
    TogglePlayback,
    ToggleGrid(bool),
    Next,
    SpeedChanged(f32),
}

impl Application for EulerFluidSimulation {
    type Message = Message;
    type Theme = Theme;
    type Executor = executor::Default;
    type Flags = ();

    fn new(_flags: ()) -> (Self, Command<Message>) {
        (
            Self {
                speed: 5,
                ..Self::default()
            },
            Command::none(),
        )
    }

    fn title(&self) -> String {
        String::from("Euler Fluid Simulation")
    }

    fn update(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::Grid(message, version) => {
                if version == self.version {
                    self.grid.update(message);
                }
            }
            Message::Tick(_) | Message::Next => {
                self.queued_ticks = (self.queued_ticks + 1).min(self.speed);

                // if let Some(task) = self.grid.tick(self.queued_ticks) {
                //     if let Some(speed) = self.next_speed.take() {
                //         self.speed = speed;
                //     }

                //     self.queued_ticks = 0;

                //     let version = self.version;

                //     return Command::perform(task, move |message| {
                //         Message::Grid(message, version)
                //     });
                // }
            }
            Message::TogglePlayback => {
                self.is_playing = !self.is_playing;
            }
            Message::ToggleGrid(show_grid_lines) => {
                self.grid.toggle_lines(show_grid_lines);
            }
            Message::SpeedChanged(speed) => {
                if self.is_playing {
                    self.next_speed = Some(speed.round() as usize);
                } else {
                    self.speed = speed.round() as usize;
                }
            }
        }

        Command::none()
    }

    fn subscription(&self) -> Subscription<Message> {
        if self.is_playing {
            time::every(Duration::from_millis(1000 / self.speed as u64))
                .map(Message::Tick)
        } else {
            Subscription::none()
        }
    }

    fn view(&self) -> Element<Message> {
        let version = self.version;
        let selected_speed = self.next_speed.unwrap_or(self.speed);
        let controls = view_controls(
            self.is_playing,
            self.grid.are_lines_visible(),
            selected_speed,
        );

        let content = column![
            self.grid
                .view()
                .map(move |message| Message::Grid(message, version)),
            controls,
        ];

        container(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    fn theme(&self) -> Theme {
        Theme::Dark
    }
}

fn view_controls<'a>(
    is_playing: bool,
    is_grid_enabled: bool,
    speed: usize,
) -> Element<'a, Message> {
    let playback_controls = row![
        button(if is_playing { "Pause" } else { "Play" })
            .on_press(Message::TogglePlayback),
        button("Next")
            .on_press(Message::Next)
            .style(theme::Button::Secondary),
    ]
    .spacing(10);

    let speed_controls = row![
        slider(1.0..=1000.0, speed as f32, Message::SpeedChanged),
        text(format!("x{speed}")).size(16),
    ]
    .width(Length::Fill)
    .align_items(Alignment::Center)
    .spacing(10);

    row![
        playback_controls,
        speed_controls,
        checkbox("Grid", is_grid_enabled, Message::ToggleGrid)
            .size(16)
            .spacing(5)
            .text_size(16)
    ]
    .padding(10)
    .spacing(20)
    .align_items(Alignment::Center)
    .into()
}