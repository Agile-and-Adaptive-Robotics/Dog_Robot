// Master Micro Testing Main Script

// This script serves to test the functionality of individual master microcontroller modules.

// Include the associated header file.
#include "Master_Micro_Header.h"

// Define global constants.
const uint8_t multiplexer_pins1[3] = {2, 3, 4};
const uint8_t multiplexer_pins2[3] = {5, 6, 7};
const uint8_t * multiplexer_port = &PORTD;

// Define global variables (local to this file). 
static volatile struct slave_struct_array slave_manager;


// Implement the main function.
int main ( void )
{
	
	// Setup the microcontroller.
	setup_micro(  );

	// Initialize the slave manager.
	initialize_slave_manager( &slave_manager );

	//// Set the multiplexer to use channel 63.
	//set_multiplexer_channel( 63 );

	// Create an empty loop.
	while( 1 ){  }

}


// Implement a function to retrieve data from Matlab, exchange command and sensory data with the slave microcontrollers, and send sensory data to Matlab at consistent time intervals.
ISR( TIMER1_COMPA_vect )						// Timer Interrupt Function 1.
{
	
	//// Read the desired pressures from Matlab.
	//usart_read_matlab_desired_pressures( &slave_manager );
	//
	//// Write the desired pressures to the slaves while collecting their most recent sensory information.
	//spi_read_write_slave_sensor_specific_command( &slave_manager, 255, 1, 255 );			// 255 = All Slaves, 1 = Command 1 (i.e., Desired Pressure), 255 = All Sensors.
	//
	//// Write the slave sensory information to Matlab.
	//usart_write_matlab_sensor_data( &slave_manager, stdout );



	//uint16_t adc_value = readADC( 0 );
	//uint16_t dac_value = ADC2DAC( adc_value );
	//write2DAC( dac_value );



	//uint16_t adc_value = readADC( 0 );
//
	//uint16_t desired_pressure_uint16 = uint102uint16( adc_value );
//
	//for ( uint8_t k = 0; k < slave_manager.num_slaves; ++k )
	//{
		//slave_manager.slave[k].desired_pressure = desired_pressure_uint16
	//}


	uint16_t desired_pressure_uint16 = 65535;

	for ( uint8_t k = 0; k < slave_manager.num_slaves; ++k )
	{
		slave_manager.slave[k].desired_pressure = desired_pressure_uint16;
	}


	// Write the desired pressures to the slaves while collecting their most recent sensory information.
	//spi_read_write_slave_sensor_specific_command( &slave_manager, 255, 1, 255 );			// 255 = All Slaves, 1 = Command 1 (i.e., Desired Pressure), 255 = All Sensors.
	spi_read_write_slave_sensor_specific_command( &slave_manager, 1, 1, 255 );			// 255 = All Slaves, 1 = Command 1 (i.e., Desired Pressure), 255 = All Sensors.


}

