#include "../JuceLibraryCode/JuceHeader.h"

#include "PerryThePlatypus.h"
#include "PerryThePlatypus_types.h"

#if JUCE_VERSION >= 0x050400
using Parameter = AudioProcessorValueTreeState::Parameter;
#endif

struct onParamChangeListener : AudioProcessorValueTreeState::Listener
{
    onParamChangeListener(PerryThePlatypusStackData* sd) : SD(sd)
    {
    }

    void parameterChanged (const String& parameterID, float newValue) override
    {
        const int idx = getParameterIndex(parameterID);
        onParamChangeCImpl(SD, idx, static_cast<double>(newValue));
    }
    
private:
    int getParameterIndex (const String& parameterID)
    {
        (void)parameterID;
    
        if (parameterID == "BYPASS") {
            return 0;
        }
        if (parameterID == "GAIN_dB") {
            return 1;
        }
        if (parameterID == "DECAY") {
            return 2;
        }
        if (parameterID == "MIX") {
            return 3;
        }
        if (parameterID == "LS_FREQ") {
            return 4;
        }
        if (parameterID == "LS_GAIN") {
            return 5;
        }
        if (parameterID == "LS_Q") {
            return 6;
        }
        if (parameterID == "HP_FREQ") {
            return 7;
        }
        if (parameterID == "HP_Q") {
            return 8;
        }
        if (parameterID == "LP_FREQ") {
            return 9;
        }
        if (parameterID == "LP_Q") {
            return 10;
        }
        if (parameterID == "LOW_PEAK") {
            return 11;
        }
        if (parameterID == "LOWP_GAIN") {
            return 12;
        }
        if (parameterID == "LOWP_Q") {
            return 13;
        }
        if (parameterID == "MID_PEAK") {
            return 14;
        }
        if (parameterID == "MIDP_GAIN") {
            return 15;
        }
        if (parameterID == "MIDP_Q") {
            return 16;
        }
        if (parameterID == "HIGH_PEAK") {
            return 17;
        }
        if (parameterID == "HIGHP_GAIN") {
            return 18;
        }
        if (parameterID == "HIGHP_Q") {
            return 19;
        }
    
        return (-1);
    }    
    
    PerryThePlatypusStackData *SD;
};

//==============================================================================
class PerryThePlatypusAudioProcessor  : public AudioProcessor
{
    //==============================================================================
#if JUCE_VERSION >= 0x050400
    const StringArray m_choices1;

public:
    PerryThePlatypusAudioProcessor()
        : paramListener(&mStackData),
          m_choices1({ "off", "on" }),
          parameters(*this, nullptr, "PerryThePlatypus", {
                std::make_unique<Parameter>("BYPASS", juce::String::fromUTF8(u8"Bypass"), juce::String::fromUTF8(u8""),
                    NormalisableRange<float>(0.f, m_choices1.size()-1.f, 1.f), 0.f,
                    [=](float value) { return m_choices1[(int) (value + 0.5)]; },
                    [=](const String& text) { return (float) m_choices1.indexOf(text); }, false, true, true),
                std::make_unique<Parameter>("GAIN_dB", juce::String::fromUTF8(u8"Gain"), juce::String::fromUTF8(u8"dB"),
NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("DECAY", juce::String::fromUTF8(u8"Decay"), juce::String::fromUTF8(u8"S"),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("MIX", juce::String::fromUTF8(u8"Mix"), juce::String::fromUTF8(u8"%"),
NormalisableRange<float>(0.f,100.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LS_FREQ", juce::String::fromUTF8(u8"LOW SHELF FREQ"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 300.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LS_GAIN", juce::String::fromUTF8(u8"LOW SHELF GAIN"), juce::String::fromUTF8(u8"dB"),
NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LS_Q", juce::String::fromUTF8(u8"LOW SHELF Q FACTOR"), juce::String::fromUTF8(u8""),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("HP_FREQ", juce::String::fromUTF8(u8"HI Pass FREQ"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 700.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("HP_Q", juce::String::fromUTF8(u8"HP Q FACTOR"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LP_FREQ", juce::String::fromUTF8(u8"Low Pass FREQ"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 14000.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LP_Q", juce::String::fromUTF8(u8"LP Q FACTOR"), juce::String::fromUTF8(u8""),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LOW_PEAK", juce::String::fromUTF8(u8"Low Freq Peak"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 50.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LOWP_GAIN", juce::String::fromUTF8(u8"Low Peak Gain"), juce::String::fromUTF8(u8"dB"),
NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("LOWP_Q", juce::String::fromUTF8(u8"LOW PEAK Q FACTOR"), juce::String::fromUTF8(u8""),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("MID_PEAK", juce::String::fromUTF8(u8"Mid Freq Peak"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 1000.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("MIDP_GAIN", juce::String::fromUTF8(u8"Mid Peak Gain"), juce::String::fromUTF8(u8"dB"),
NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("MIDP_Q", juce::String::fromUTF8(u8"MID Peak Q FACTOR"), juce::String::fromUTF8(u8""),
NormalisableRange<float>(0.5f,2.f), 0.5f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("HIGH_PEAK", juce::String::fromUTF8(u8"Peak Freq"), juce::String::fromUTF8(u8"Hz"),
NormalisableRange<float>(20.f,20000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 5000.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("HIGHP_GAIN", juce::String::fromUTF8(u8"Peak Gain"), juce::String::fromUTF8(u8"dB"),
NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("HIGHP_Q", juce::String::fromUTF8(u8"Peak Q FACTOR"), juce::String::fromUTF8(u8""),
NormalisableRange<float>(0.f,1.f), 0.5f, [](float val) {return String(val, 3);}, nullptr) })

    {
        mStackData.pd = &mPersistentData;
        
        PerryThePlatypus_initialize(&mStackData);

        createPluginInstance(&mStackData, reinterpret_cast<unsigned long long>(this));

        parameters.addParameterListener("BYPASS", &paramListener);
        parameters.addParameterListener("GAIN_dB", &paramListener);
        parameters.addParameterListener("DECAY", &paramListener);
        parameters.addParameterListener("MIX", &paramListener);
        parameters.addParameterListener("LS_FREQ", &paramListener);
        parameters.addParameterListener("LS_GAIN", &paramListener);
        parameters.addParameterListener("LS_Q", &paramListener);
        parameters.addParameterListener("HP_FREQ", &paramListener);
        parameters.addParameterListener("HP_Q", &paramListener);
        parameters.addParameterListener("LP_FREQ", &paramListener);
        parameters.addParameterListener("LP_Q", &paramListener);
        parameters.addParameterListener("LOW_PEAK", &paramListener);
        parameters.addParameterListener("LOWP_GAIN", &paramListener);
        parameters.addParameterListener("LOWP_Q", &paramListener);
        parameters.addParameterListener("MID_PEAK", &paramListener);
        parameters.addParameterListener("MIDP_GAIN", &paramListener);
        parameters.addParameterListener("MIDP_Q", &paramListener);
        parameters.addParameterListener("HIGH_PEAK", &paramListener);
        parameters.addParameterListener("HIGHP_GAIN", &paramListener);
        parameters.addParameterListener("HIGHP_Q", &paramListener);

    }
    //==============================================================================
#else // For JUCE prior to 5.4.0
public:
    PerryThePlatypusAudioProcessor()
    :   paramListener(&mStackData), parameters (*this, nullptr)
    {
        mStackData.pd = &mPersistentData;
        
        PerryThePlatypus_initialize(&mStackData);

        createPluginInstance(&mStackData, reinterpret_cast<unsigned long long>(this));

        //
        // Parameter property BYPASS
        //
        const StringArray choices1({ "off", "on" });
        parameters.createAndAddParameter ("BYPASS", juce::String::fromUTF8(u8"Bypass"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.f, choices1.size()-1.f, 1.f), 0.f,
            [=](float value) { return choices1[(int) (value + 0.5)]; },
            [=](const String& text) { return (float) choices1.indexOf(text); },
            false, true, true);
        parameters.addParameterListener("BYPASS", &paramListener);

        //
        // Parameter property GAIN_dB
        //
        parameters.createAndAddParameter ("GAIN_dB", juce::String::fromUTF8(u8"Gain"), juce::String::fromUTF8(u8"dB"),
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("GAIN_dB", &paramListener);

        //
        // Parameter property DECAY
        //
        parameters.createAndAddParameter ("DECAY", juce::String::fromUTF8(u8"Decay"), juce::String::fromUTF8(u8"S"),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("DECAY", &paramListener);

        //
        // Parameter property MIX
        //
        parameters.createAndAddParameter ("MIX", juce::String::fromUTF8(u8"Mix"), juce::String::fromUTF8(u8"%"),
            NormalisableRange<float>(0.f, 100.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("MIX", &paramListener);

        //
        // Parameter property LS_FREQ
        //
        parameters.createAndAddParameter ("LS_FREQ", juce::String::fromUTF8(u8"LOW SHELF FREQ"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            300.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LS_FREQ", &paramListener);

        //
        // Parameter property LS_GAIN
        //
        parameters.createAndAddParameter ("LS_GAIN", juce::String::fromUTF8(u8"LOW SHELF GAIN"), juce::String::fromUTF8(u8"dB"),
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LS_GAIN", &paramListener);

        //
        // Parameter property LS_Q
        //
        parameters.createAndAddParameter ("LS_Q", juce::String::fromUTF8(u8"LOW SHELF Q FACTOR"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LS_Q", &paramListener);

        //
        // Parameter property HP_FREQ
        //
        parameters.createAndAddParameter ("HP_FREQ", juce::String::fromUTF8(u8"HI Pass FREQ"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            700.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("HP_FREQ", &paramListener);

        //
        // Parameter property HP_Q
        //
        parameters.createAndAddParameter ("HP_Q", juce::String::fromUTF8(u8"HP Q FACTOR"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("HP_Q", &paramListener);

        //
        // Parameter property LP_FREQ
        //
        parameters.createAndAddParameter ("LP_FREQ", juce::String::fromUTF8(u8"Low Pass FREQ"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            14000.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LP_FREQ", &paramListener);

        //
        // Parameter property LP_Q
        //
        parameters.createAndAddParameter ("LP_Q", juce::String::fromUTF8(u8"LP Q FACTOR"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LP_Q", &paramListener);

        //
        // Parameter property LOW_PEAK
        //
        parameters.createAndAddParameter ("LOW_PEAK", juce::String::fromUTF8(u8"Low Freq Peak"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            50.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LOW_PEAK", &paramListener);

        //
        // Parameter property LOWP_GAIN
        //
        parameters.createAndAddParameter ("LOWP_GAIN", juce::String::fromUTF8(u8"Low Peak Gain"), juce::String::fromUTF8(u8"dB"),
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LOWP_GAIN", &paramListener);

        //
        // Parameter property LOWP_Q
        //
        parameters.createAndAddParameter ("LOWP_Q", juce::String::fromUTF8(u8"LOW PEAK Q FACTOR"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("LOWP_Q", &paramListener);

        //
        // Parameter property MID_PEAK
        //
        parameters.createAndAddParameter ("MID_PEAK", juce::String::fromUTF8(u8"Mid Freq Peak"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            1000.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("MID_PEAK", &paramListener);

        //
        // Parameter property MIDP_GAIN
        //
        parameters.createAndAddParameter ("MIDP_GAIN", juce::String::fromUTF8(u8"Mid Peak Gain"), juce::String::fromUTF8(u8"dB"),
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("MIDP_GAIN", &paramListener);

        //
        // Parameter property MIDP_Q
        //
        parameters.createAndAddParameter ("MIDP_Q", juce::String::fromUTF8(u8"MID Peak Q FACTOR"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.5f, 2.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("MIDP_Q", &paramListener);

        //
        // Parameter property HIGH_PEAK
        //
        parameters.createAndAddParameter ("HIGH_PEAK", juce::String::fromUTF8(u8"Peak Freq"), juce::String::fromUTF8(u8"Hz"),
            NormalisableRange<float>(20.f, 20000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            5000.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("HIGH_PEAK", &paramListener);

        //
        // Parameter property HIGHP_GAIN
        //
        parameters.createAndAddParameter ("HIGHP_GAIN", juce::String::fromUTF8(u8"Peak Gain"), juce::String::fromUTF8(u8"dB"),
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("HIGHP_GAIN", &paramListener);

        //
        // Parameter property HIGHP_Q
        //
        parameters.createAndAddParameter ("HIGHP_Q", juce::String::fromUTF8(u8"Peak Q FACTOR"), juce::String::fromUTF8(u8""),
            NormalisableRange<float>(0.f, 1.f), 0.5f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("HIGHP_Q", &paramListener);

        parameters.state = ValueTree(Identifier("PerryThePlatypus"));
    }
#endif

    //==============================================================================
    ~PerryThePlatypusAudioProcessor()
    {
        PerryThePlatypus_terminate(&mStackData);
    }
    
    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override
    {
        (void)samplesPerBlock;
        resetCImpl(&mStackData, sampleRate);
        setLatencySamples(getLatencyInSamplesCImpl(&mStackData));
    }

    void releaseResources() override                { }
    
    
    void processBlock (AudioBuffer<double>& buffer, MidiBuffer&) override
    {
        ScopedNoDenormals noDenormals;
        const int nSamples = buffer.getNumSamples();
        const int nChannels = buffer.getNumChannels();
        const double* const* pin = buffer.getArrayOfReadPointers();
        double* const* pout = buffer.getArrayOfWritePointers();
        std::vector<const double*> inputs(pin,pin+nChannels);
        std::vector<double*> outputs(pout,pout+nChannels);
        int i_;

        PerryThePlatypusStackData *SD = &mStackData;

        if (nChannels < 2) {
            const int numExtraChannels = 2 - nChannels;
            tempBuffer.setSize(numExtraChannels, nSamples);
            if (nChannels < 2) {
                tempBuffer.clear(0, nSamples);
                const double* const* p = tempBuffer.getArrayOfReadPointers();
                std::copy(p, p+numExtraChannels, std::back_inserter(inputs));
            }
            if (nChannels < 2) {
                double* const* p = tempBuffer.getArrayOfWritePointers();
                std::copy(p, p+numExtraChannels, std::back_inserter(outputs));
            }
        }

        int osz0_;
        int osz1_;
        if (nSamples <= MAX_SAMPLES_PER_FRAME) {
            /* Fast path for common frame sizes. */
            const int isz0_ = nSamples;
            const int isz1_ = nSamples;
            processEntryPoint(SD, (double)nSamples,
                    inputs[0], &isz0_,
                    inputs[1], &isz1_,
                    outputs[0], &osz0_,
                    outputs[1], &osz1_);
        } else {
            /* Fallback for unusually large frames. */
            int isz0_ = MAX_SAMPLES_PER_FRAME;
            int isz1_ = MAX_SAMPLES_PER_FRAME;
            int n = MAX_SAMPLES_PER_FRAME;
            for (i_ = 0; i_ < nSamples; i_ += MAX_SAMPLES_PER_FRAME) {
                if (i_ + MAX_SAMPLES_PER_FRAME > nSamples) {
                    n = nSamples - i_;
                    isz0_ = nSamples - i_;
                    isz1_ = nSamples - i_;
                }
                processEntryPoint(SD, (double)n,
                        inputs[0]+i_, &isz0_,
                        inputs[1]+i_, &isz1_,
                        outputs[0]+i_, &osz0_,
                        outputs[1]+i_, &osz1_);
            }
        }

    }
    
    void processBlock (AudioBuffer<float>& buffer,  MidiBuffer& midiMessages) override
    {
        AudioBuffer<double> doubleBuffer;
        doubleBuffer.makeCopyOf(buffer);
        processBlock(doubleBuffer, midiMessages);
        buffer.makeCopyOf(doubleBuffer);
    }
    
    //==============================================================================
    bool hasEditor() const override                 { return true; }
    AudioProcessorEditor* createEditor() override;
    
    //==============================================================================
    const String getName() const override           { return JucePlugin_Name; }

    bool acceptsMidi() const override               { return false; }
    bool producesMidi() const override              { return false; }
    bool isMidiEffect () const override             { return false; }
    double getTailLengthSeconds() const override    { return 0.0;   }

    //==============================================================================
    // NB: some hosts don't cope very well if you tell them there are 0 programs,
    // so this should be at least 1, even if you're not really implementing programs.
    int getNumPrograms() override                       { return 1;  }
    int getCurrentProgram() override                    { return 0;  }
    void setCurrentProgram (int index) override         { (void) index; }
    const String getProgramName (int index) override    { (void) index; return {}; }
    void changeProgramName (int index, const String& newName) override  { (void) index; (void) newName; }
    
    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override
    {
        auto xml (parameters.state.createXml());
        copyXmlToBinary (*xml, destData);
    }
    
    void setStateInformation (const void* data, int sizeInBytes) override
    {
        auto xmlState (getXmlFromBinary (data, sizeInBytes));
        if (xmlState != nullptr)
            if (xmlState->hasTagName (parameters.state.getType()))
                parameters.state = ValueTree::fromXml (*xmlState);
    }
    
    bool supportsDoublePrecisionProcessing() const override  { return true; }
    
private:
    //==============================================================================
    static const int MAX_SAMPLES_PER_FRAME = 4096;

    PerryThePlatypusStackData mStackData;
    PerryThePlatypusPersistentData mPersistentData;
    onParamChangeListener paramListener;
    AudioBuffer<double> tempBuffer;
    
    //==============================================================================
    AudioProcessorValueTreeState parameters;
 
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PerryThePlatypusAudioProcessor)
};

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PerryThePlatypusAudioProcessor();
}

#include "PerryThePlatypusPluginEditor.h"

AudioProcessorEditor* PerryThePlatypusAudioProcessor::createEditor()
{
    return new PerryThePlatypusAudioProcessorEditor(*this, parameters);
}

