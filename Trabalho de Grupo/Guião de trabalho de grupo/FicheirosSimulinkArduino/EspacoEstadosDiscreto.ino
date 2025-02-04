//Definição de nº de estados e entradas

#define STATES 6  // Number of states
#define INPUTS 2  // Number of inputs
//Definição da Matriz A

float A[STATES][STATES] = {
{0.043996,-0.52853,0.44611,-0.0027668,0.005226,0.012053},
{-0.52862,0.49,-0.081764,0.0052273,0.0092865,0.017263},
{0.44624,-0.081584,-0.036882,0.01205,0.017265,0.014504},
{47.227,-5.1154,-30.035,0.044151,-0.52831,0.4463},
{-5.1419,17.197,-34.975,-0.52902,0.48943,-0.08225},
{-29.943,-35.017,12.11,0.44765,-0.079562,-0.035183}
};

// MatrizB -

float B[STATES][INPUTS] = {
{0.0056194,0.0012005},
{0.0099855,0.0017778},
{0.018565,0.0012758},
{-0.56808,0.057656},
{-0.549,0.00066875},
{-0.085551,-0.11286}
};


// Kr

float Kr[STATES] =  {16.624,110.3,-55.716,0.51505,-0.43561,0.10403};

// initial values for states
float x[STATES] =  {0.0, 0.0, 0.0, 0.0, 0.0,0.0};
                 
// Initial state
float x_next[STATES] =  {0.0, 0.0, 0.0, 0.0, 0.0,0.0};
float u[INPUTS] =  {0.0, 0.0};

#define DT 50 // em milisegundos
unsigned long tActual = 0; unsigned long tInit = 0;

float uControl = 0.0;
float yAtual = 0.0;

unsigned long tIni;
unsigned long t;
void setup()
{
  // put your setup code here, to run once:
  Serial.begin(115200);
  pinMode(11, OUTPUT); //PWM PIN 11  with PWM wire
  TCCR2B = TCCR2B & B11111000 | B00000001; // for PWM frequency of 31372.55 Hz
   tIni = millis();
}



void loop() {
 
  tInit = millis();
   // read the value from the sensor:
  int sensorValue = analogRead(A5);
 
     updateState(); // update state
    
     uControl = computeControl(); // compute the next system input
    
     while (millis() <= tInit + 2) {
     // wait for next sample
      }

}




void matrixMultiply(float *result, float matrix[][STATES], float *vector, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

void matrixMultiplyOut(float *result, float matrix[][INPUTS], float *vector, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}
void updateState() {
    tInit = millis();
    // Compute A * x
    float Ax[STATES];
    matrixMultiply(Ax, A, x, STATES, STATES);

    // Compute B * u
    float Bu[STATES];
    matrixMultiplyOut(Bu, B, u, STATES, INPUTS);

    // Update x_next = Ax + Bu
    for (int i = 0; i < STATES; i++) {
        x_next[i] = Ax[i] + Bu[i];
    }

    // Update the state
    for (int i = 0; i < STATES; i++) {
        x[i] = x_next[i];
    }
    
}


float computeControl() {
    // Compute -Kr * x
    float Control=0;
    for (int i = 0; i < STATES; i++) {
        Control += Kr[i]*x[i];
    }
    /*
    for (int i = 0; i < STATES; i++) {
       Serial.print(x[i]);
       Serial.print(" ");
   }
    
    for (int i = 0; i < STATES; i++) {
       Serial.print(Kr[i]);
       Serial.print(" ");
   }
     Serial.println(Control);
     */
    return Control;
    
}
