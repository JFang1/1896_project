#include <Adafruit_Sensor.h>

#include <Wire.h>
#include <SD.h>
#include <Adafruit_BNO055.h>

/* This driver uses the Adafruit unified sensor library (Adafruit_Sensor),
   which provides a common 'type' for sensor data and some helper functions.
   To use this driver you will also need to download the Adafruit_Sensor
   library and include it in your libraries folder.
   You should also assign a unique ID to this sensor for use with
   the Adafruit Sensor API so that you can identify this particular
   sensor in any data logs, etc.  To assign a unique ID, simply
   provide an appropriate value in the constructor below (12345
   is used by default in this example).
   Connections
   ===========
   Connect SCL to analog 5
   Connect SDA to analog 4
   Connect VDD to 3-5V DC
   Connect GROUND to common ground
   History
   =======
   2015/MAR/03  - First release (KTOWN)
   2015/AUG/27  - Added calibration and system status helpers
*/

/* Set the delay between fresh samples */
#define BNO055_SAMPLERATE_DELAY_MS (10)

Adafruit_BNO055 bno = Adafruit_BNO055(55);
String dataString = "";
File dataFile;
const int led = 8;
const int button = 3;
uint8_t sys = 0;
uint8_t gyro = 0;
uint8_t accel = 0;
uint8_t mag=0;


/**************************************************************************/
/*
    Display some basic info about the sensor status
*/
/**************************************************************************/
void displaySensorStatus(void)
{
  /* Get the system status values (mostly for debugging purposes) */
  uint8_t system_status, self_test_results, system_error;
  system_status = self_test_results = system_error = 0;
  bno.getSystemStatus(&system_status, &self_test_results, &system_error);

  /* Display the results in the Serial Monitor */
  Serial.println("");
  Serial.print("System Status: 0x");
  Serial.println(system_status, HEX);
  Serial.print("Self Test:     0x");
  Serial.println(self_test_results, HEX);
  Serial.print("System Error:  0x");
  Serial.println(system_error, HEX);
  Serial.println("");
  delay(500);
}

/**************************************************************************/
/*
    Display sensor calibration status
*/
/**************************************************************************/
void displayCalStatus(void)
{
  /* Get the four calibration values (0..3) */
  /* Any sensor data reporting 0 should be ignored, */
  /* 3 means 'fully calibrated" */
  sys = gyro = accel = mag = 0;
  bno.getCalibration(&sys, &gyro, &accel, &mag);

  // Display the individual values if connected to PC
  Serial.print("Sys:");
  Serial.print(sys, DEC);
  Serial.print(" G:");
  Serial.print(gyro, DEC);
  Serial.print(" A:");
  Serial.print(accel, DEC);
  Serial.print(" M:");
  Serial.print(mag, DEC);

  //When the sensor is fully calibrated, make the LED flash
  if(mag==3 && gyro ==3 && accel==3){
    delay(200);
    digitalWrite(led,HIGH);
    delay(200);
    digitalWrite(led,LOW);
  }

 //Print new line before next calibration output
 Serial.print("\r\n");
  
}

/**************************************************************************/
/*
    Arduino setup function (automatically called at startup)
*/
/**************************************************************************/
void setup(void)
{
  Serial.begin(9600);
  Serial.println("Orientation Sensor Test"); Serial.println("");

  /* Initialise the sensor */
  if(!bno.begin())
  {
    /* There was a problem detecting the BNO055 ... check your connections */
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    while(1);
  }

  /* Optional: Display current status */
  displaySensorStatus();

  /* Set operation mode to NDOF */
  Adafruit_BNO055::adafruit_bno055_opmode_t mode;
  mode = Adafruit_BNO055::adafruit_bno055_opmode_t::OPERATION_MODE_NDOF;
  bno.setMode(mode);
  bno.setExtCrystalUse(true);


   // Open serial communications and wait for port to open:
  Serial.begin(9600);
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB port only
  }
  
  Serial.print("Initializing SD card...");

  //Alerts us if the initialization of the SD card fails
  if (!SD.begin(10)) {
    Serial.println("initialization failed!");
    return;
  }

  //Alert when SD card successfully initialized
  Serial.println("initialization done.");
  
  
  SD.remove("testcal.txt");

  //Open File for data output
  dataFile = SD.open("testcal.txt", FILE_WRITE);

  //Set the pinModes for the led and button
  pinMode(led,OUTPUT);
  pinMode(button,INPUT);

  //Infinite loop used to calibrate sensor.
  //Can break out of loop by pressing the button set above
  while(true){
    displayCalStatus();
    if(digitalRead(button)==HIGH){
      break;  
    }
  }

  //Turn the LED off to acknowledge calibration loop has been exited
  digitalWrite(led,LOW);
  
  //Wait 4 seconds before starting data collection
  delay(4000);

  //Turn LED on to signal data collection has begun
  digitalWrite(led,HIGH);
 
}

/**************************************************************************/
/*
    Arduino loop function, called once 'setup' is complete (your own code
    should go here)
*/
/**************************************************************************/
void loop(void)
{
  //Delay 9ms to guarantee that we are not logging data at a rate greater than 100Hz
  delay(9);
  
  // Get a new sensor event
  sensors_event_t event;
  bno.getEventAcc(&event);

  //Get calibration data to add to data string
  bno.getCalibration(&sys, &gyro, &accel, &mag);

  //Build Data String for Output File
  dataString+=event.acceleration.x;
  dataString+=",";
  dataString+=event.acceleration.y;
  dataString+=",";
  dataString+=event.acceleration.z;
  dataString+=",";
  dataString+=millis();
  dataString+=",";
  dataString+=(sys);
  dataString+="\r\n";

  //Write Data String to File
  dataFile.println(dataString);

  //Reset Data String to null before building next sample
  dataString = ""; 

  //When button is pressed, turn off LED and close output file
  //Once entered, we cannot exit this loop
  if(digitalRead(button)==HIGH){
    while(true){
      if(dataFile){
        dataFile.close();
      }
      digitalWrite(led,LOW);
    }
  }
  
}
