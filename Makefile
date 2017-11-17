CC = mpiifort
CFLAGS =
TARGET = LJ
OBJECTS = basics.o math_module.o mc_module.o LJ.f90

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

math_module.o : math_module.f90
	$(CC) $(CFLAGS) -c $^

basics.o : basics.f90
	$(CC) $(CFLAGS) -c $^

mc_module.o : mc_module.f90 basics.o
	$(CC) $(CFLAGS) -I basics.o -c mc_module.f90

clean :
	rm *.o *.mod *.dat LJ
