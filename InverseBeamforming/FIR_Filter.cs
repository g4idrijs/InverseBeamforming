using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public class FIR_Filter
	{
		private double[] _filter;

		public FIR_Filter(double[] filter)
		{
			this._filter = filter;
		}

		/// <summary>
		/// Implements FIR filtering on the signal using the coefficients given when the instance was constructed
		/// </summary>
		/// <param name="signal">Signal to filter</param>
		/// <param name="filt">Use the RF or DS filter (RF by default)</param>
		/// <returns>Filtered signal</returns>
		public double[] Filter(double[] signal, double[] filter=null)
		{
			int numTaps = _filter.Length;
			int numSigPoints = signal.Length;
			int top = 0, n, k;

			double[] reg = new double[numTaps];
			double[] filteredSignal = new double[numSigPoints];
			double y;

			for (int j = 0; j < numSigPoints; j++)
			{
				reg[top] = signal[j];
				y = 0;
				n = 0;

				for (k = top; k >= 0; k--)
				{
					y += _filter[n++] * reg[k];
				}
				for (k = numTaps - 1; k > top; k--)
				{
					y += _filter[n++] * reg[k];
				}
				filteredSignal[j] = y;
				top++;
				if (top >= numTaps)
				{
					top = 0;
				}
			}

			return filteredSignal;
		}

		public void set_taps(double[] taps)
		{
			_filter = taps;
		}
	} 

}
