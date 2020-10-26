#!/bin/bash

TOP=`pwd`

rm -fr ${TOP}/out.*

cd ${TOP}/stm32_cubemx/nucleo_l476rg
make clean

cd ${TOP}
