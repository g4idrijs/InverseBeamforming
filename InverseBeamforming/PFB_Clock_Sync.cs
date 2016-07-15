using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public class PFB_Clock_Sync
	{
		private bool _updated;
		private double _sps;
		private double _sample_num;
		private double _loop_bw;
		private double _damping;
		private double _alpha;
		private double _beta;

		private int _nfilters;
		private int _taps_per_filter;
		private FIR_Filter[] _filters;
		private FIR_Filter[] _diff_filters;
		private List<List<double>> _taps;
		private List<List<double>> _dtaps;
		private double[] _update_taps;
		Complex f;

		private double _k;
		private double _rate;
		private double _rate_i;
		private double _rate_f;
		private double _max_dev;
		private int _filtnum;
		private int _osps;
		private double _error;
		private int _out_idx;

		private UInt64 _ol_in, _new_in, _last_out;

		/// <summary>
		/// Construct a new instance of the PFB Clock Sync class
		/// </summary>
		/// <param name="sps">Samples per symbol. Note: This is not the same samples per symbol used in the Modulation classes. Once I figure out what this number is, I will update this</param> //TODO: Update what sps means in Gnuradio
		/// <param name="loop_bw">Loop bandwidth of the filters. This parameter sets alpha and beta</param>
		/// <param name="taps">Taps for each of the filters</param>
		/// <param name="filter_size">Number of filters. The default should be 32</param>
		/// <param name="init_phase">Initial starting phase</param>
		/// <param name="max_rate_deviation">How far the rate can deviate from zero</param>
		/// <param name="osps">Number of output samples per symbol</param>
		public PFB_Clock_Sync(double sps, double loop_bw, double[] taps, int filter_size, double init_phase, double max_rate_deviation, int osps)
		{
			_updated = false;
			_nfilters = filter_size;
			_max_dev = max_rate_deviation;
			_osps = osps;
			_error = 0;
			_out_idx = 0;

			if (taps.Length == 0)
				throw new ArgumentException("PFB_Clock_Sync: please specify a filter.");

			_nfilters = filter_size;
			_sps = Math.Floor(sps);

			// Set the damping factor for a critically damped system
			_damping = 2 * _nfilters;

			// Set the bandwidth, which will then call update_gains()
			set_loop_bandwidth(loop_bw);

			// Store the last filter between calls to work
			// The accumulator keeps track of overflow to increment the stride correctly.
			// set it here to the fractional difference based on the initial phaes
			_k = init_phase;
			_rate = (sps - Math.Floor(sps)) * (double)_nfilters;
			_rate_i = (int)Math.Floor(_rate);
			_rate_f = _rate - (float)_rate_i;
			_filtnum = (int)Math.Floor(_k);

			_filters = new FIR_Filter[_nfilters];
			_diff_filters = new FIR_Filter[_nfilters];

			// Create an FIR filter for each channel and zero out the taps
			double[] vtaps = null;
			for (int i = 0; i < _nfilters; i++)
			{
				_filters[i] = new FIR_Filter(vtaps);
				_diff_filters[i] = new FIR_Filter(vtaps);
			}

			// Now, actually set the filters' taps
			List<double> dtaps = null;
			create_diff_taps(taps, out dtaps);
			set_taps(taps, out _taps, ref _filters);
			set_taps(dtaps.ToArray(), out _dtaps, ref _diff_filters);

			_ol_in = 0;
			_new_in = 0;
			_last_out = 0;
		}

		//TODO: Finish implementing this PFB clock sync code (it could be tough, there might be a need to come up with a replacement for the tag system)
		public int general_work(int noutput_items, Complex[] input, ref Complex[] output)
		{
			if (_updated)
			{
				List<double> dtaps;
				create_diff_taps(_update_taps, out dtaps);
				set_taps(_update_taps, out _taps, ref _filters);
				set_taps(dtaps.ToArray(), out _dtaps, ref _diff_filters);
				_updated = false;
				return 0;            // history requirements may have changed.
			}

			double[] err = null, outrate = null, outk = null;
			if (input.Length == 4)
			{
				err = (double*)output_items[1];
				outrate = (double*)output_items[2];
				outk = (double*)output_items[3];
			}

			int i = 0, count = 0;
			double error_r, error_i;

			// produce output as long as we can and there are enough input samples
			while (i < noutput_items)
			{
				if (tags.size() > 0)
				{
					size_t offset = tags[0].offset - nitems_read(0);
					if ((offset >= (size_t)count) && (offset < (size_t)(count + _sps)))
					{
						double center = (double)pmt::to_double(tags[0].value);
						_k = _nfilters * (center + (offset - count));

						tags.erase(tags.begin());
					}
				}

				while (_out_idx < _osps)
				{

					_filtnum = (int)Math.Floor(_k);

					// Keep the current filter number in [0, _nfilters]
					// If we've run beyond the last filter, wrap around and go to next sample
					// If we've gone below 0, wrap around and go to previous sample
					while (_filtnum >= _nfilters)
					{
						_k -= _nfilters;
						_filtnum -= _nfilters;
						count += 1;
					}
					while (_filtnum < 0)
					{
						_k += _nfilters;
						_filtnum += _nfilters;
						count -= 1;
					}

					output[i + _out_idx] = _filters[_filtnum].Filter(input[count + _out_idx]);
					_k = _k + _rate_i + _rate_f; // update phase


					// Manage Tags
					std::vector<tag_t> xtags;
					std::vector<tag_t>::iterator itags;
					_new_in = nitems_read(0) + count + _out_idx + _sps;

					get_tags_in_range(xtags, 0, _ol_in, _new_in);
					for (itags = xtags.begin(); itags != xtags.end(); itags++)
					{
						tag_t new_tag = *itags;
						//new_tag.offset = _last_out + _taps_per_filter/(2*_sps) - 2;
						new_tag.offset = _last_out + _taps_per_filter / 4 - 2;

						ad_item_tag(0, new_tag);
					}
					_ol_in = _new_in;
					_last_out = nitems_written(0) + i + _out_idx;

					_out_idx++;

					if (output_items.size() == 4)
					{
						err[i] = _error;
						outrate[i] = _rate_f;
						outk[i] = _k;
					}

					// We've run out of output items we can create; return now.
					if (i + _out_idx >= noutput_items)
					{
						return i;
					}
				}

				// reset here; if we didn't complete a full osps samples last time,
				// the early return would take care of it.
				_out_idx = 0;

				// Update the phase and rate estimates for this symbol
				Complex diff = _diff_filters[_filtnum].Filter(input[count]);
				error_r = output[i].Real * diff.Real;
				error_i = output[i].Imaginary * diff.Imaginary;
				_error = (error_i + error_r) / 2.0;       // average error from I&Q channel

				// Run the control loop to update the current phase (k) and
				// tracking rate estimates based on the error value
				// Interpolating here to update rates for ever sps.
				for (int s = 0; s < _sps; s++)
				{
					_rate_f = _rate_f + _beta * _error;
					_k = _k + _rate_f + _alpha * _error;
				}

				// Keep our rate within a good range
				_rate_f = branchlessClip(_rate_f, _max_dev);

				i += _osps;
				count += (int)Math.Floor(_sps);
			}
			
			return i;
		}






		/// <summary>
		/// Set the loop bandwidth and update alpha and beta
		/// </summary>
		/// <param name="bw">New loop bandwidth</param>
		void set_loop_bandwidth(double bw)
		{
			if (bw < 0)
			{
				throw new ArgumentOutOfRangeException("pfb_clock_sync_ccf: invalid bandwidth. Must be >= 0.");
			}

			_loop_bw = bw;
			update_gains();
		}

		/// <summary>
		/// Update alpha and beta
		/// </summary>
		void update_gains()
		{
			double denom = (1.0 + 2.0 * _damping * _loop_bw + _loop_bw * _loop_bw);
			_alpha = (4 * _damping * _loop_bw) / denom;
			_beta = (4 * _loop_bw * _loop_bw) / denom;
		}

		/// <summary>
		/// Get the taps needed for the filtering
		/// </summary>
		/// <param name="newtaps"></param>
		/// <param name="ourtaps"></param>
		/// <param name="ourfilter"></param>
		void set_taps(double[] newtaps, out List<List<double>> ourtaps, ref FIR_Filter[] ourfilter)
		{
			int i, j;

			int ntaps = newtaps.Length;
			_taps_per_filter = (int)Math.Ceiling((double)ntaps / (double)_nfilters);

			// Create _numchan vectors to store each channel's taps
			ourtaps = new List<List<double>>();
			for (int k = 0; k < _nfilters; k++)
			{
				ourtaps.Add(new List<double>());
			}

			// Make a vector of the taps plus fill it out with 0's to fill
			// each polyphase filter with exactly _taps_per_filter
			List<double> tmp_taps;
			tmp_taps = newtaps.ToList();
			while ((float)(tmp_taps.Count()) < _nfilters * _taps_per_filter)
			{
				tmp_taps.Add(0);
			}

			// Partition the filter
			for (i = 0; i < _nfilters; i++)
			{
				// Each channel uses all _taps_per_filter with 0's if not enough taps to fill out
				for (j = 0; j < _taps_per_filter; j++)
				{
					ourtaps[i].Add(tmp_taps[i + j * _nfilters]);
				}

				// Build a filter for each channel and add it's taps to it
				ourfilter[i].set_taps(ourtaps[i].ToArray());
			}
		}

		/// <summary>
		/// Create the differential filter taps
		/// </summary>
		/// <param name="newtaps">taps used to create the differential filter</param>
		/// <param name="difftaps">Taps for the differential filter</param>
		void create_diff_taps(double[] newtaps, out List<double> difftaps)
		{
			double[] diff_filter = new double[3];
			diff_filter[0] = -1;
			diff_filter[1] = 0;
			diff_filter[2] = 1;

			double pwr = 0;
			double tap = 0;
			difftaps = new List<double>();
			difftaps.Add(0);
			for (int i = 0; i < newtaps.Length - 2; i++)
			{
				tap = 0;
				for (int j = 0; j < diff_filter.Length; j++)
				{
					tap += diff_filter[j] * newtaps[i + j];
				}
				difftaps.Add(tap);
				pwr += Math.Abs(tap);
			}
			difftaps.Add(0);

			// Normalize the taps
			for (int i = 0; i < difftaps.Count; i++)
			{
				difftaps[i] *= _nfilters / pwr;
				if (difftaps[i] != difftaps[i])
				{
					throw new InvalidOperationException("pfb_clock_sync_ccf::create_diff_taps produced NaN.");
				}
			}
		}
	}
}
