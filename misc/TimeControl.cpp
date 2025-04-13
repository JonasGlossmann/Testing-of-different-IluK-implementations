//
// Created by jonas on 27.11.2024.
//

#include "TimeControl.h"

#include <iostream>

TimeControl::TimeControl() {
    this->startTime = std::chrono::high_resolution_clock::now();
    line = 0;
}

void TimeControl::checkout() {
    this->endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = this->endTime - this->startTime;
    std::cout <<"Number Timecheck: "<< this->line << " Time needed: " << elapsed.count() << " ms" << std::endl;
    this->line++;
    this->startTime = std::chrono::high_resolution_clock::now();
}

