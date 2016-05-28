using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public static class ExtensionMethods
	{
		/// <summary>
		/// Returns the mean value of the array
		/// </summary>
		/// <param name="doubleArray">Array to find the mean of</param>
		/// <returns>Mean value of all elements in the array</returns>
		public static double Mean(this double[] doubleArray)
		{
			double mean = 0;

			for(int i=0; i<doubleArray.Length; i++)
			{
				mean += doubleArray[i];
			}

			return mean / doubleArray.Length;
		}

		/// <summary>
		/// Returns the variance of the values in the array
		/// </summary>
		/// <param name="doubleArray">Array to find the variance of</param>
		/// <returns>Variance of the values in the array</returns>
		public static double Variance(this double[] doubleArray)
		{
			double variance = 0;
			double mean = doubleArray.Mean();

			for (int i = 0; i < doubleArray.Length; i++)
			{
				variance += Math.Pow(doubleArray[i] - mean, 2);
			}

			return variance / doubleArray.Length;
		}

		/// <summary>
		/// Returns the variance of the values in the array
		/// </summary>
		/// <param name="doubleArray">Array to find the variance of</param>
		/// <param name="mean">Mean of the array values</param>
		/// <returns>Variance of the values in the array</returns>
		public static double Variance(this double[] doubleArray, double mean)
		{
			double variance = 0;

			for (int i = 0; i < doubleArray.Length; i++)
			{
				variance += Math.Pow(doubleArray[i] - mean, 2);
			}

			return variance / doubleArray.Length;
		}

		/// <summary>
		/// Returns the standard deviation of the values in the array
		/// </summary>
		/// <param name="doubleArray">Array to find the standard deviation of</param>
		/// <returns>Standard deviation of the values in the array</returns>
		public static double StandardDeviation(this double[] doubleArray)
		{
			//The standard deviation is the square root of the variance
			return Math.Sqrt(doubleArray.Variance());
		}

		/// <summary>
		/// Returns the standard deviation of the array when the mean is already known
		/// </summary>
		/// <param name="doubleArray">Array to find the standard deviation of</param>
		/// <param name="mean">Mean value of the elements in the array</param>
		/// <returns>Standard deviation of the elements in the array</returns>
		public static double StandardDeviation(this double[] doubleArray, double mean)
		{
			//The standard deviation is just the square root of the variance
			return Math.Sqrt(doubleArray.Variance(mean));
		}
	}
}
