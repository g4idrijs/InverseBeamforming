using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// Contains all the components necesary to perform simulations regaurding digital communications systems
/// and inverse beamforming
/// </summary>
namespace InverseBeamforming
{
	/// <summary>
	/// Provides useful extension methods for common use cases
	/// </summary>
	public static class ExtensionMethods
	{
		private static readonly double CPLXDMATH_ZERO_TEST = 1.0E-50;

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

		/// <summary>
		/// Clones the array assuming it isn't a null reference
		/// </summary>
		/// <typeparam name="T">Type of the array to clone</typeparam>
		/// <param name="arr">Array to clone</param>
		/// <returns>Null if arr is a null reference, or a clone of the array if it is not a null reference</returns>
		public static T[] NullCloneSafely<T>(this T[] arr)
		{
			if(arr==null)
			{
				return null;
			}
			else
			{
				return (T[])arr.Clone();
			}
		}

		/// <summary>
		/// Calls the OnNext method of each simulation reporter in the list
		/// </summary>
		/// <param name="list">List of simulation reporters</param>
		/// <param name="sim">Simulation to report on</param>
		public static void OnNext_List(this List<Simulations.SimulationReporter> list, Simulations.SingleSimulation sim)
		{
			foreach (var item in list)
			{
				item.OnNext(sim);
			}
		}

		/// <summary>
		/// Calls the OnStart method of each simulation reporter in the list
		/// </summary>
		/// <param name="list">List of simulation reporters</param>
		public static void OnStart_List(this List<Simulations.SimulationReporter> list)
		{
			foreach (var item in list)
			{
				item.OnStart();
			}
		}

		/// <summary>
		/// Calls the OnError method of each simulation reporter in the list
		/// </summary>
		/// <param name="list">List of simulation reporters</param>
		/// <param name="e">Exception that caused the error</param>
		public static void OnError_List(this List<Simulations.SimulationReporter> list, Exception e)
		{
			foreach (var item in list)
			{
				item.OnError(e);
			}
		}

		/// <summary>
		/// Calls the OnCompleted method of each simulation reporter in the list
		/// </summary>
		/// <param name="list">List of simulation reporters</param>
		public static void OnCompleted_List(this List<Simulations.SimulationReporter> list)
		{
			foreach (var item in list)
			{
				item.OnCompleted();
			}
		}

		/// <summary>
		/// Calls the Subscribe method of each simulation report in the list for the given provider
		/// </summary>
		/// <param name="list">List of simulation reporters</param>
		/// <param name="provider">Simulation tracker for the simulation</param>
		public static void Subscribe_List(this List<Simulations.SimulationReporter> list, Simulations.SimulationTracker provider)
		{
			foreach (var item in list)
			{
				item.Subscribe(provider);
			}
		}

		/// <summary>
		/// Clears the contents of the array
		/// </summary>
		/// <param name="arr">Array to clear</param>
		public static void Clear(this double[] arr)
		{
			for(int i=0; i<arr.Length; i++)
			{
				arr[i]= 0;
			}
		}

		/// <summary>
		/// Copies startValue into each element of the array
		/// </summary>
		/// <param name="arr">Array to fill with startValue</param>
		/// <param name="startValue">Value to put in each element of the array</param>
		public static void Fill(this double[] arr, double startValue = 0)
		{
			for (int i = 0; i < arr.Length; i++)
			{
				arr[i] = startValue;
			}
		}

		/// <summary>
		/// Take the square root of the double
		/// </summary>
		/// <param name="val">Double to take the square root of</param>
		/// <returns>Square root of the given value</returns>
		public static double Sqrt(this double val)
		{
			return Math.Sqrt(val);
		}

		/// <summary>
		/// Raise a complex number to a power
		/// </summary>
		/// <param name="X">Complex number</param>
		/// <param name="B">Power</param>
		/// <returns>Complex number to that power</returns>
		public static Complex Pow(this Complex X, double B)
		{
			double Theta, Radius;
			if (Math.Abs(B) < CPLXDMATH_ZERO_TEST) // anything, inc 0, to the 0th pow = 1
			{
				return new Complex(1.0, 0.0);
			}
			if (X.Magnitude < CPLXDMATH_ZERO_TEST) // 0 to any power, except 0, = 0
			{
				return new Complex(0.0, 0.0);
			}
			Theta = B * Math.Atan2(X.Imaginary, X.Real); // atan2(Y,X) Realturns +/- PI/2 if X = 0  errs with X=Y=0
			Radius = Math.Pow(Math.Sqrt(X.Real * X.Real + X.Imaginary * X.Imaginary), B);
			return new Complex(Radius * Math.Cos(Theta), Radius * Math.Sin(Theta));
		}

		/// <summary>
		/// Cosine of a Complex number
		/// </summary>
		/// <param name="X">Complex nuber to take the cosine of</param>
		/// <returns>Cosine of the number</returns>
		public static Complex Cos(this Complex X)
		{
			return new Complex(Math.Cos(X.Real) * Math.Cosh(X.Imaginary), -Math.Sin(X.Real) * Math.Sinh(X.Imaginary));
		}

		/// <summary>
		/// Returns the conjugate of the complex number
		/// </summary>
		/// <param name="X">Complex number oto take the conjugate of</param>
		/// <returns>Conjugate of the complex number</returns>
		public static Complex Conjugate(this Complex X)
		{
			return new Complex(X.Real, -X.Imaginary);
		}
	}
}
