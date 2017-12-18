main.exe: main.cpp
	g++ -I/usr/include/eigen3/ main.cpp -o main.exe -lsfml-graphics -lsfml-system -lsfml-window
