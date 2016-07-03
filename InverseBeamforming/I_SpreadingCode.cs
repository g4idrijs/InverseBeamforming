using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// Contains all the components necesary to perform simulations regaurding digital communications systems
/// and inverse beamforming
/// </summary>
namespace InverseBeamforming
{
	/// <summary>
	/// Collects many different waveform modulations into 1 class
	/// </summary>
	public partial class Modulations
	{
		/// <summary>
		/// Defines an interface to control the way spreading codes are generated and handled.
		/// </summary>
		interface I_SpreadingCode
		{
			/// <summary>
			/// Holds the spreading code matrix
			/// </summary>
			byte[,] CodeMatrix { get; set; }

			/// <summary>
			/// Number of chips per symbol
			/// </summary>
			int NumChips { get; set; }

			/// <summary>
			/// Gets a full waveform length spreading code, ready to multiply with the signal. 
			/// </summary>
			/// <param name="user">Index of the spreading code to use (zero based)</param>
			/// <param name="numSamples">Number of samples per symbol duration</param>
			/// <param name="numChips">Number of chips per symbol duration</param>
			/// <param name="numSymbols">Number of symbols</param>
			/// <returns>Array containing the spreading code ready to be mixed with the signal</returns>
			byte[] GetSpreadingCode(int user, int numSymbols);

			/// <summary>
			/// The number of samples in one symbol period (Samples/symbol)
			/// </summary>
			int SamplesPerSymbol { get; set; }

			/// <summary>
			/// Number of symbols in each waveform
			/// </summary>
			int NumberSymbolsPerWaveform { get; set; }

			/// <summary>
			/// Applies a spreading code to a waveform. The original waveform is modified
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="numSamples">Number of samples per symbol</param>
			/// <param name="numChips">Number of chips per symbol</param>
			/// <param name="numSymbols">Number of symbols in the waveform</param>
			void SpreadWaveform(ref double[] waveform, int user);

			/// <summary>
			/// Applies a spreading code to a waveform. The original waveform is modified
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="numSamples">Number of samples per symbol</param>
			/// <param name="numChips">Number of chips per symbol</param>
			/// <param name="numSymbols">Number of symbols in the waveform</param>
			/// <remarks>In the implementation, it is functionally equivalent to SpreadWaveform. Applying the same spreading code twice gives the original signal back.</remarks>
			void DespreadWaveform(ref double[] waveform, int user);
		}
	}
}
