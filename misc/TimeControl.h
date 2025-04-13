//
// Created by jonas on 27.11.2024.
//

#ifndef TIMECONTROL_H
#define TIMECONTROL_H
#include <chrono>


class TimeControl {
    int line;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime ;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime ;
    public:
    TimeControl();
    void checkout();

};



#endif //TIMECONTROL_H
