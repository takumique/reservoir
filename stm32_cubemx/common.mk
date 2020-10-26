C_SOURCES += \
$(wildcard $(APP_PATH)/*.c) \
$(wildcard $(APP_PATH)/cmcis/*.c)

C_DEFS += \
-DPRECISION_F32 \
-DCONST_WEIGHTS

C_INCLUDES += \
-I$(APP_PATH) \
-I$(APP_PATH)/cmcis

LDFLAGS += \
-u _printf_float

LIBS += \
-larm_cortexM4lf_math

LIBDIR += \
-LDrivers/CMSIS/DSP/Lib/GCC
