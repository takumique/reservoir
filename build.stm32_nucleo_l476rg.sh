#!/bin/bash

TOP=`pwd`

export APP_PATH=${TOP}/app
export GCC_PATH=/usr/bin

cd ${TOP}/stm32_cubemx/nucleo_l476rg
make

cd ${TOP}
