To compile:
``g++ -o main read.cpp 4401-81.cpp adch.cpp mina.cpp main.cpp`` 

To run:
``./main ``

In main.cpp fill bomb parametrs, initial motion parameters, parameters of corrective engine, target position(if need)

To get the affected area:
Fill all parameters without target position, compile and run.

To get parametrs of trajectory at specific target:
Fill all parameters. Uncomment ``cout`` in function ``update()``, ``get_res()`` in ``otrMathModel.cpp``. Comment all after first `` res = mina.get_res(dt, iter); `` in ``main.cpp``. Compile and run.   