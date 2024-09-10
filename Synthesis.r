require 'open-uri'
require 'mspire'
require 'mspire/mzml'
require 'csv'

if __FILE__ == $0


  # do not include any '/' or '_' in any sample names
  files = Dir["C:/Users/cindy/Downloads/lab 2021-2022/synthetic-data-ruby/*_intensity*"]

  for i in 0..files.length-1
      files[i] = files[i].gsub(/.*\//, '').partition("_int")[0]
  end

  for x in 0..files.length-1
      rttable = CSV.read(files[x] + "_rt.csv")
      mztable = CSV.read(files[x] + "_mz.csv")
      inttable = CSV.read(files[x] + "_intensity.csv")

      spects = []
      for i in 0..rttable[0].length-1

        time = rttable[0][i]
        mz = []
        int =[]
        for j in 0..mztable[i].length-1
          mz << mztable[i][j].to_f
          int << inttable[i][j].to_f
        end

        spec1 = Mspire::Mzml::Spectrum.new('scan=1') do |spec|
          # profile and ms_level 1
          # http://repo.msys2.org/distrib/x86_64/msys2-x86_64-20190524.exe
          spec.describe_many!(['MS:1000128', ['MS:1000511', 1]])
          spec.data_arrays = [
            Mspire::Mzml::DataArray.new(mz).describe!('MS:1000514'), # mz
            Mspire::Mzml::DataArray.new(int).describe!('MS:1000515')  # intensity
          ]
          spec.scan_list = Mspire::Mzml::ScanList.new do |sl|
            scan = Mspire::Mzml::Scan.new do |scan|
              # retention time of 42 seconds
              scan.describe! 'MS:1000016', time, 'UO:0000010'     # rt
            end
            sl << scan
          end
        end

        spects << spec1
      end

      mzml = Mspire::Mzml.new do |mzml|
        mzml.id = 'ms1'
        mzml.cvs = Mspire::Mzml::CV::DEFAULT_CVS
        mzml.file_description = Mspire::Mzml::FileDescription.new  do |fd|
          fd.file_content = Mspire::Mzml::FileContent.new
          fd.source_files << Mspire::Mzml::SourceFile.new
        end
        default_instrument_config = Mspire::Mzml::InstrumentConfiguration.new("IC").describe!('MS:1000031')
        mzml.instrument_configurations << default_instrument_config
        software = Mspire::Mzml::Software.new
        mzml.software_list << software
        default_data_processing = Mspire::Mzml::DataProcessing.new("id")
        mzml.data_processing_list << default_data_processing
        mzml.run = Mspire::Mzml::Run.new("run", default_instrument_config) do |run|
          spectrum_list = Mspire::Mzml::SpectrumList.new(default_data_processing, spects) # [spec1, spec2]
          run.spectrum_list = spectrum_list
        end
      end

      file_type = ".mzML"
      mzml.write(files[x] + file_type)
      # mg.say_bye
  end

end

# writes ms2 spectras: 
# spec2 = Mspire::Mzml::Spectrum.new('scan=2') do |spec|
#     # centroid,  ms_level 2, MSn spectrum,
#     spec.describe_many!(['MS:1000127', ['MS:1000511', 2], "MS:1000580"])
#     spec.data_arrays = [
#       Mspire::Mzml::DataArray[1,2,3.5].describe!('MS:1000514'),
#       Mspire::Mzml::DataArray[5,6,5].describe!('MS:1000515')
#     ]
#     spec.scan_list = Mspire::Mzml::ScanList.new do |sl|
#       scan = Mspire::Mzml::Scan.new do |scan|
#         # retention time of 42 seconds
#         scan.describe! 'MS:1000016', 45.0, 'UO:0000010'
#       end
#       sl << scan
#     end
#     precursor = Mspire::Mzml::Precursor.new( spec1.id )
#     si = Mspire::Mzml::SelectedIon.new
#     # the selected ion m/z:
#     si.describe! "MS:1000744", 2.0
#     # the selected ion charge state
#     si.describe! "MS:1000041", 2
#     # the selected ion intensity
#     si.describe! "MS:1000042", 5
#     precursor.selected_ions = [si]
#     spec.precursors = [precursor]
#   end