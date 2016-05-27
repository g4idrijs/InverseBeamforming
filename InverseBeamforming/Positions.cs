using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Provides methods to read positions, calculate distances and time delays based on the speed of light.
	/// </summary>
    public class Positions
    {
		/// <summary>
		/// Positions of transmit antennas, coordinates given in meters
		/// </summary>
		public double[,] TxPositions
		{
			get
			{
				return this._txPositions;
			}
			set
			{
				this._txPositions = value;
			}
		}
		private double[,] _txPositions;

		/// <summary>
		/// Positions of the Receive antennas, coordinates given in meters
		/// </summary>
		public double[,] RxPositions
		{
			get
			{
				return this._rxPositions;
			}
			set
			{
				this._rxPositions = value;
			}
		}
		private double[,] _rxPositions;

		/// <summary>
		/// Constructs a new Positions class with the given TX and RX positions
		/// </summary>
		/// <param name="txPositions">2D array of the positions of the transmit antennas</param>
		/// <param name="rxPositions">2D array of the positions of the receive antennas</param>
		public Positions(double[,] txPositions, double[,] rxPositions)
		{
			this._txPositions = txPositions;
			this._rxPositions = rxPositions;
		}

		/// <summary>
		/// Finds the distance between one receive antenna, 0 indexed and specified by rxIndex, to every transmit antenna.
		/// </summary>
		/// <param name="rxIndex">Zero based index of the receive antenna of interest</param>
		/// <returns>Array of distances from the desired receive antenna to each of the transmit antennas
		/// NOTE: Can return null</returns>
		public double[] getDistances(int rxIndex)
		{
			//There are no positions
			if (this._rxPositions == null || this._txPositions == null)
				return null;

			//The dimensions of the rx and tx positions do not match
			if (this._txPositions.GetUpperBound(1) != this._rxPositions.GetUpperBound(1))
				return null;
			
			//Get the number of dimensions in the transmitter positions
			var numRows = this._txPositions.GetUpperBound(0) + 1;

			//Get the total number of transmitters
			var numDims = this._txPositions.GetUpperBound(1) + 1;

			//Initialize and array of distances
			var distances = new double[numRows];

			
			double sum;	

			//Loop through each transmit antenna
			for (int i = 0; i < numRows; i++)
			{
				//Reset sum
				sum = 0;

				//Sum the squared distances in each dimension
				for(int k=0; k < numDims; k++)
				{
					sum += Math.Pow(_rxPositions[rxIndex, k] - _txPositions[i, k],2);
				}
				//Take the square root to find the overall distance
				distances[i] = Math.Sqrt(sum);
			}

			//Return the distance
			return distances;
		}

		/// <summary>
		/// Returns each the distance between each pair of transmit antennas and receive antennas.
		/// </summary>
		/// <returns>Distances between each pair of transmit and receive antenna. 
		/// Format: Distances in columns come from same receive antenna, rows come from same transmit antenna.
		/// NOTE: Can return null</returns>
		public double[,] getAllDistances()
		{
			//There are no positions
			if (this._rxPositions == null || this._txPositions == null)
				return null;

			//The dimensions of the rx and tx positions do not match
			if (this._txPositions.GetUpperBound(1) != this._rxPositions.GetUpperBound(1))
				return null;

			//Get the number of Transmit and receive antennas
			var numTx = this._txPositions.GetUpperBound(0) + 1;
			var numRx = this._rxPositions.GetUpperBound(0) + 1;

			//Get the total number of transmitters
			var numDims = this._txPositions.GetUpperBound(1) + 1;

			//Intialize the array
			double[,] distances = new double[numTx, numRx];

			double sum;

			//Loop through each receive antenna
			for(int j=0; j<numRx; j++)
			{
				//Loop through each transmit antenna
				for (int i = 0; i < numTx; i++)
				{
					//Reset sum
					sum = 0;

					//Sum the squared distances in each dimension
					for (int k = 0; k < numDims; k++)
					{
						sum += Math.Pow(_rxPositions[j, k] - _txPositions[i, k], 2);
					}
					//Take the square root to find the overall distance
					distances[i,j] = Math.Sqrt(sum);
				}
			}

			//Return the array of distances
			return distances;
		}

		/// <summary>
		/// Calculates the time delays on signals between each transmit antenna and a receive antenna of interest
		/// </summary>
		/// <param name="rxIndex">Receive antenna of interest</param>
		/// <returns>Returns an array of time delays representing the amount of time required for the signal to propagate from each TX antenna to the receive antenna of interest
		/// NOTE: Can return null</returns>
		public double[] getTimeDelays(int rxIndex)
		{
			//Get the distances from each transmitter to the receiver of interest
			double[] distances = getDistances(rxIndex);

			//If the distances function didn't return anything, then return null as well;
			if (distances == null)
				return null;

			//Define an array to hold the time delays. Should be the same length as the distances array
			double[] timeDelays = new double[distances.Length];

			//Define the speed of light
			const double c = 299792458;

			//Loop through each distance
			for (int i=0; i<distances.Length; i++)
			{
				//Time delay is the distance traveled over the speed of light
				timeDelays[i] = distances[i] / c;
			}

			//Return the time delays
			return timeDelays;
		}

		/// <summary>
		/// Calculate the time delay of signals between each pair of a transmit and receive antenna
		/// </summary>
		/// <returns>2D array of all the time delays between each pair of transmit and receive antennas.
		/// Format: Time delays in columns come from same receive antenna, rows come from same transmit antenna.
		/// NOTE: Returns null if unsuccessful.</returns>
		public double[,] getAllTimeDelays()
		{
			//Get the distances from each transmitter to the receiver of interest
			double[,] distances = getAllDistances();

			//If the distances function didn't return anything, then return null as well;
			if (distances == null)
				return null;

			//Get the number of Transmit and receive antennas
			var numTx = this._txPositions.GetUpperBound(0) + 1;
			var numRx = this._rxPositions.GetUpperBound(0) + 1;

			//Define an array to hold the time delays. Should be the same length as the distances array
			double[,] timeDelays = new double[numTx, numRx];

			//Define the speed of light
			const double c = 299792458;

			//Loop through each distance
			for (int i = 0; i < numTx; i++)
			{
				for(int k=0; k< numRx; k++)
				{
					//Time delay is the distance traveled over the speed of light
					timeDelays[i,k] = distances[i,k] / c;
				} 
			}

			//Return the time delays
			return timeDelays;
		}

		/// <summary>
		/// Calculates the sample delays to place on each link in order to compensate for different distances between antennas
		/// </summary>
		/// <param name="rxIndex">Receiver antenna of interest</param>
		/// <param name="samplingRate">Sampling rate of the signal</param>
		/// <returns>An array containing the number of samples to delay each signal by according to the distance between the pair.
		/// NOTE: Returns null if not successful.</returns>
		public int[] getSampleDelays(int rxIndex, double samplingRate)
		{
			//Get the time delays
			double[] timeDelays = getTimeDelays(rxIndex);

			//If the time delays are null, then return null here as well
			if (timeDelays == null)
				return null;

			//Initialize an array to hold the sample delays
			int[] sampleDelays = new int[timeDelays.Length];

			//Loop through each pair and calculate the sample delay
			for(int i=0; i<timeDelays.Length; i++)
			{
				sampleDelays[i] = (int)(timeDelays[i] / samplingRate);
			}

			//Return the sample delays
			return sampleDelays;
		}

		/// <summary>
		/// Calculates the sample delay for each pair of transmit and receive antennas
		/// </summary>
		/// <param name="samplingRate">Sampling rate of the signal</param>
		/// <returns>2D array containing the number of samples to delay each signal by according to the distance between the pair.
		/// Format: Time delays in columns come from same receive antenna, rows come from same transmit antenna.
		/// NOTE: Returns null if not successful.</returns>
		public int[,] getAllSampleDelays(double samplingRate)
		{
			//Get the time delays
			double[,] timeDelays = getAllTimeDelays();

			//If the time delays are null, then return null here as well
			if (timeDelays == null)
				return null;

			//Get the number of Transmit and receive antennas
			var numTx = this._txPositions.GetUpperBound(0) + 1;
			var numRx = this._rxPositions.GetUpperBound(0) + 1;

			//Initialize an array to hold the sample delays
			int[,] sampleDelays = new int[numTx,numRx];

			//Loop through each pair and calculate the sample delay
			for (int i = 0; i < numTx; i++)
			{
				//Loop through each receive antenna
				for (int k = 0; k < numRx; k++)
				{
					sampleDelays[i, k] = (int)(timeDelays[i, k] / samplingRate);
				} 
			}

			//Return the sample delays
			return sampleDelays;
		}
	}
}
