all:
	g++ -Ofast optymalizacja.cpp -o optymalizacja -std=gnu++17 -lboost_random -lboost_system \
	-Wall -Wextra -Wpedantic -Wshadow -Wenum-compare -Wunreachable-code -Werror=narrowing -Werror=return-type
run:
	optymalizacja
