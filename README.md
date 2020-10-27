# reservoir
This project is to translate [this example](https://qiita.com/pokotsun/items/dd8eb48fadeee052110b) to C and run on tiny MCU.  
## Getting started
### Prerequisite packages
```
sudo apt install -y cmake
```
### Linux machine
Run following command to build:
```
./build.sh
```
Run following command to run:
```
./out.generic/app/reservoir
```
Make sure training.txt is in your current directory if you undeclare CONST_WEIGHTS in app/CMakeLists.txt .
predict.txt is generated in your current directory.

### STM32 MCU
Run following command to build:
```
./build.sh -t <target device>
```
Currently only supported target device is "stm32_nucleo_l476rg" . The build artifacts are generated under stm32_cubemx/<target device>/build directory.  
Run following command to flash:
```
st-flash --format ihex write ./stm32_cubemx/nucleo_l476rg/build/nucleo_l476rg.hex
```
To evaluate target device is running correctly, use UART2 via ST-LINK 2-1 to see logs. In case your ST-LINK virtual serial console is on /dev/ttyACM0 , run following command to see logs.
```
minicom -D /dev/ttyACM0
```
You will see following logs and green LED on target board is light if everything went well.
```
Welcome to minicom 2.7.1

OPTIONS: I18n 
Compiled on Dec 23 2019, 02:06:26.
Port /dev/ttyACM0, 08:45:49

Press CTRL-A Z for help on special keys

INIT
res->in_weights: -0.100000 0.100000 0.100000 0.100000
res->res_weights: 0.000245 -0.013844 0.027784 -0.028285
res->out_weights: 0.000000 0.000000 0.000000 0.000000
res->res_nodes: 0.000000 0.000000 0.000000 0.000000
res->x: 0.000000 0.000000 0.000000 0.000000
res->y: 0.000000 0.000000 0.000000 0.000000
TRAINED
res->in_weights: -0.100000 0.100000 0.100000 0.100000
res->res_weights: 0.000245 -0.013844 0.027784 -0.028285
res->out_weights: 0.177775 0.102471 -0.224610 0.327216
res->res_nodes: 0.038026 -0.033332 -0.038418 -0.041508
res->x: 0.000000 0.000000 0.000000 0.000000
res->y: 0.000000 0.000000 0.000000 0.000000
```
## Tips
### Installing ST-Link
```
sudo apt install -y git cmake libusb-1.0
(configure git)
git clone https://github.com/texane/stlink.git
make
cd build/Release
sudo make install
```
Generally libs are installed under /usr/local/lib . Therefore, make sure your /etc/ld.conf.d/ is correctly search libs under /usr/local/lib .  
```
sudo ldconfig
```
to reload libs.
### Installing minicom
```
sudo apt install -y minicom
```


