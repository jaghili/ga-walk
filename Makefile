EIGENDIR = /usr/include/eigen3/
LFLAGS = -lsfml-graphics -lsfml-system -lsfml-window

main.exe: main.cpp
	g++ -I$(EIGENDIR) $< -o $@ $(LFLAGS)

test.exe: test.cpp
	g++ -I$(EIGENDIR) $< -o $@ $(LFLAGS)
