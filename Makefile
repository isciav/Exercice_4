CC = g++
CFLAGS = -Wall -pedantic -ansi -g 
EXEC_NAME = Exercice4
INCLUDES =
LIBS =
OBJ_FILES = Exercice4.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES) *.out

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

Exercice4.o: ConfigFile.tpp
