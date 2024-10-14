classdef PerryThePlatypus < audioPlugin
    % Parametric EQ with Reverb
    % Aditya Jeganath K
  
    properties
        
        BYPASS = 'off';

        GAIN_dB = 0; %dB

        %EQ Parameters 
        
        %High Shelf
        HS_FREQ = 10000;     
        HS_GAIN = 0; 
        HS_Q = 0.5; 

        LS_FREQ = 300;  % Hz 
        LS_GAIN = 0; %dB
        LS_Q = 0.5;
        
        LP_FREQ = 14000;
        LP_Q = 0.5;

        HP_FREQ = 700;% Hz 
        HP_Q = 0.5; 
        
        LOW_PEAK = 50;% Hz 
        LOWP_GAIN = 0; %dB
        LOWP_Q = 0.5;
        
        MID_PEAK = 1000;% Hz 
        MIDP_GAIN = 0; %dB
        MIDP_Q = 0.5;

        HIGH_PEAK = 5000;% Hz 
        HIGHP_GAIN = 0; %dB
        HIGHP_Q = 0.5;
        
        %Reverb Parameters
        DECAY = 0.5;
        MIX = 0;
    end
%% Plugin Interface
    properties (Constant)

        PluginInterface = audioPluginInterface( ...
             'PluginName','PerrythePlatypus',...
             'VendorName','TropicThunder',...
             'VendorVersion','1.0.0',...
             'UniqueId','TopG',...
          audioPluginGridLayout(...
            'RowHeight', [40 100 25 100 40 100],...
            'ColumnWidth', [80 10 80 10 80 40 80 10 80 40 80 40 80 40 80 ]),...  
          audioPluginParameter('BYPASS',...
            'DisplayName', 'Bypass',...
            'Mapping', {'enum', 'off', 'on'},...
            'Layout', [2,15],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'switch_metal.png',...
            'FilmstripFrameSize', [64 64]),...
          audioPluginParameter('GAIN_dB',...
            'Label', 'dB',...
            'DisplayName', 'Gain',...
            'Mapping',{'lin',-12,12},'Style', 'rotaryknob',...
            'Layout', [4,15],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),... 
          audioPluginParameter('DECAY',...
            'DisplayName','Decay',...
            'Label','S',...
            'Mapping',{'lin',0,1},'Style', 'rotaryknob',...
            'Layout',[2,13],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('MIX',...
            'DisplayName','Mix',...
            'Label','%',...
            'Mapping',{'lin',0,100},'Style', 'rotaryknob',...
            'Layout', [4,13],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LS_FREQ',...
            'DisplayName','LOW SHELF FREQ',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,11],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LS_GAIN',...
            'DisplayName','LOW SHELF GAIN',...
            'Label', 'dB',...
            'Mapping',{'lin',-12,12},'Style', 'rotaryknob',...
            'Layout', [4,11],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LS_Q',...
            'DisplayName','LOW SHELF Q FACTOR',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [6,11],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
           audioPluginParameter('HP_FREQ',...
            'DisplayName','HI Pass FREQ',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,9],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('HP_Q',...
            'DisplayName','HP Q FACTOR',...
            'Label', 'Hz',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [4,9],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LP_FREQ',...
            'DisplayName','Low Pass FREQ',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,7],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LP_Q',...
            'DisplayName','LP Q FACTOR',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [4,7],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LOW_PEAK',...
            'DisplayName','Low Freq Peak',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,5],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LOWP_GAIN',...
            'DisplayName','Low Peak Gain',...
            'Label', 'dB',...
            'Mapping',{'lin',-12,12},'Style', 'rotaryknob',...
            'Layout', [4,5],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('LOWP_Q',...
            'DisplayName','LOW PEAK Q FACTOR',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [6,5],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('MID_PEAK',...
            'DisplayName','Mid Freq Peak',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,3],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
          audioPluginParameter('MIDP_GAIN',...
            'DisplayName','Mid Peak Gain',...
            'Label', 'dB',...
            'Mapping',{'lin',-12,12},'Style', 'rotaryknob',...
            'Layout', [4,3],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
           audioPluginParameter('MIDP_Q',...
            'DisplayName','MID Peak Q FACTOR',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [6,3],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...  
           audioPluginParameter('HIGH_PEAK',...
            'DisplayName','Peak Freq',...
            'Label', 'Hz',...
            'Mapping',{'log',20,20000},'Style', 'rotaryknob',...
            'Layout', [2,1],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
           audioPluginParameter('HIGHP_GAIN',...
            'DisplayName','Peak Gain',...
            'Label', 'dB',...
            'Mapping',{'lin',-12,12},'Style', 'rotaryknob',...
            'Layout', [4,1],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100]),...
           audioPluginParameter('HIGHP_Q',...
            'DisplayName','Peak Q FACTOR',...
            'Mapping',{'lin',0.5,2},'Style', 'rotaryknob',...
            'Layout', [6,1],...
            'DisplayNameLocation', 'above',...
            'Filmstrip', 'knob.png',...
            'FilmstripFrameSize', [100  100])...
        );
    end

     properties (Access = private)
      
        %filter_HS = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_LS = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);     
        filter_LP = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_HP = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_HPF = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_MPF = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_LPF = struct('w', zeros(2, 2), 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);

        % parameter smoothers for frequency ranges 
        
        hsFreqSmoother
        lsFreqSmoother
        hpFreqSmoother
        lpFreqSmoother
        lpfSmoother
        mpfSmoother
        hpfSmoother

        % parameter smoothers for gain 
        hsGainSmoother
        lsGainSmoother
        
        lpGainSmoother
        lpfGainSmoother
        mpfGainSmoother
        hpfGainSmoother

        % parameter smoothers for Q factor 
        hsQSmoother
        lsQSmoother
        hpQSmoother
        lpQSmoother
        lpfQSmoother
        mpfQSmoother
        hpfQSmoother

        GAIN_dBSmoother

        fs = 44100;
        fn=22050;
        reverb;

        gainLinear;
     end
%% Initialization
     methods
         function plugin = PerryThePlatypus()

            % Initialize frequency smoothers with default values 
            %plugin.hsFreqSmoother = ParameterSmoother(plugin.HS_FREQ, 0.9);
            plugin.lsFreqSmoother = ParameterSmoother(plugin.LS_FREQ,0.2);
            plugin.lpFreqSmoother = ParameterSmoother(plugin.LP_FREQ,0.2);
            plugin.hpFreqSmoother = ParameterSmoother(plugin.HP_FREQ,0.2);
            plugin.lpfSmoother = ParameterSmoother(plugin.LOW_PEAK,0.2);
            plugin.mpfSmoother = ParameterSmoother(plugin.MID_PEAK,0.2);
            plugin.hpfSmoother = ParameterSmoother(plugin.HIGH_PEAK,0.2);

            % Initialize gain smoothers with default values
            plugin.hsGainSmoother = ParameterSmoother(plugin.HS_GAIN, 0.2);
            plugin.lsGainSmoother = ParameterSmoother(plugin.LS_GAIN, 0.2);
            plugin.lpfGainSmoother = ParameterSmoother(plugin.LOWP_GAIN, 0.2);
            plugin.mpfGainSmoother = ParameterSmoother(plugin.MIDP_GAIN, 0.2);
            plugin.hpfGainSmoother = ParameterSmoother(plugin.HIGHP_GAIN, 0.2);

            %Initialize Q smoothers with default values 
            %plugin.hsQSmoother = ParameterSmoother(plugin.HS_Q, 0.9);
            plugin.lsQSmoother = ParameterSmoother(plugin.LS_Q, 0.9);
            plugin.lpQSmoother = ParameterSmoother(plugin.LP_Q, 0.9);
            plugin.hpQSmoother = ParameterSmoother(plugin.HP_Q, 0.9);
            plugin.hpfQSmoother = ParameterSmoother(plugin.HIGHP_Q, 0.9);
            plugin.lpfQSmoother = ParameterSmoother(plugin.LOWP_Q, 0.9);
            plugin.mpfQSmoother = ParameterSmoother(plugin.MIDP_Q, 0.9);

            %Initialize Gain smoothers with default values
            plugin.GAIN_dBSmoother = ParameterSmoother(plugin.GAIN_dB, 0.9);

            %Initialize Reverb
            plugin.reverb = reverberator();
            updateReverb(plugin);

            %Initialize EQ
            
            %update_HS(plugin);
            %update_LP(plugin);
            %update_HP(plugin);
            %update_LS(plugin);
            %update_HPF(plugin);
            %update_MPF(plugin);
            %update_LPF(plugin);

            %Initialize Gain
            updateGain(plugin);
        
        end

        function out = process(plugin, in)

            out = coder.nullcopy(zeros(size(in)));
            out1 = coder.nullcopy(zeros(size(in)));
            out2 = coder.nullcopy(zeros(size(in)));

            % Smooth Value
           %updateGain(plugin);
           %update_HS(plugin);
           %update_LP(plugin);
           %update_HP(plugin);
           %update_LS(plugin);
           %update_HPF(plugin);
           %update_MPF(plugin); 
           %update_LPF(plugin);

            numChannels = min(size(in));

            for ch = 1:numChannels
                x = in(:,ch); % Input Buffer for this channel

                % Process each filter
                [x, plugin.filter_LP.w(:,ch)]  = processEQ(plugin,x, plugin.filter_LP, ch);
                [x, plugin.filter_LS.w(:,ch)]  = processEQ(plugin,x, plugin.filter_LS, ch);  
                [x, plugin.filter_HP.w(:,ch)]  = processEQ(plugin,x, plugin.filter_HP, ch);
                [x, plugin.filter_HPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_HPF, ch);
                [x, plugin.filter_LPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_LPF, ch);
                [x, plugin.filter_MPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_MPF, ch);
%{
                [x_LP, plugin.filter_LP.w(:,ch)]  = processEQ(plugin,x, plugin.filter_LP, ch);
                [x_LS, plugin.filter_LS.w(:,ch)]  = processEQ(plugin,x, plugin.filter_LS, ch);  
                [x_HP, plugin.filter_HP.w(:,ch)]  = processEQ(plugin,x, plugin.filter_HP, ch);
                [x_HPF, plugin.filter_HPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_HPF, ch);
                %[x, plugin.filter_LPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_LPF, ch);
                %[x, plugin.filter_MPF.w(:,ch)] = processEQ(plugin,x, plugin.filter_MPF, ch);
                out1(:,ch) = x_LP./4 + x_LS./4 + x_HP./4 + x_HPF./4; 
%}
                out1(:,ch) = x;
            end
            
            out1(:,:) = plugin.gainLinear * out1(:,:);          

            out2(:,:) = plugin.reverb(out1);

            if strcmp(plugin.BYPASS, 'on')
                out = in(:,:);
            else
                out(:,:) = out2(:,:);
            end

        end

             % Nested function for processing EQ
    function [y, w] = processEQ(plugin,x, filter, ch)
    y = zeros(size(x));  % Always initialize y
    w = 0;
    w1 = filter.w(1, ch); % w(n-1)
    w2 = filter.w(2, ch); 

    b0 = filter.b0;

    if abs(b0) < 1e-6
       b0 = 1e-6; 
    end
    
    a0 = filter.a0 / b0;
    a1 = filter.a1 / b0;
    a2 = filter.a2 / b0;
    b1 = filter.b1 / b0;
    b2 = filter.b2 / b0;

    for n = 1:length(x)
        w = x(n) - b1 * w1 - b2 * w2;
        y(n) = a0 * w + a1 * w1 + a2 * w2;
        
        w2 = w1;
        w1 = w;
    end

    filter.w(:,ch) = w1;
    filter.w(:,ch) = w2;
end

 %% Plugin Rest

        function reset(plugin)
            plugin.fs = getSampleRate(plugin);

            plugin.fn = plugin.fs/2;

            numChannels = size(plugin.filter_LP.w,2);
            
            %plugin.filter_HS.w = [0 0; 0 0];

            plugin.filter_LS.w = zeros(2, numChannels);

            plugin.filter_LP.w = zeros(2, numChannels);

            plugin.filter_HPF.w = zeros(2, numChannels);

            plugin.filter_HP.w = zeros(2, numChannels);

            %plugin.filter_MPF.w = [0 0; 0 0];

            %plugin.filter_LPF.w = [0 0; 0 0];
        end

        %% High Shelf Update & Set Function
%{
        function update_HS(plugin)
        f0 = plugin.HS_FREQ; %hsFreqSmoother.step();
        Q = plugin.HS_Q; %hsQSmoother.step();

        gain = plugin.hsGainSmoother.step();%hsGainSmoother.step();
        w0=2*pi*f0/plugin.fs;
        alpha=sin(w0)/(2*Q);
        A=sqrt(db2mag(gain));
            
        plugin.filter_HS.a0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha );
        plugin.filter_HS.a1 = -2*A*( (A-1) - (A+1)*cos(w0)                   );
        plugin.filter_HS.a2 =    A*( (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha );
        plugin.filter_HS.b0 =        (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha;
        plugin.filter_HS.b1 =    2*( (A-1) - (A+1)*cos(w0)                   );
        plugin.filter_HS.b2 =        (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha;
        end

        function set.HS_FREQ(plugin,val)
        plugin.HS_FREQ = val;
        plugin.hsFreqSmoother.setTargetValue(val);     %#ok   
        update_HS(plugin);
        end

        function set.HS_GAIN(plugin,val)
        plugin.HS_GAIN = val;
        plugin.hsGainSmoother.setTargetValue(val);    %#ok
        update_HS(plugin);
        end

        function set.HS_Q(plugin,val)
        plugin.HS_Q = val;    
        plugin.hsQSmoother.setTargetValue(val);
        update_HS(plugin);
        end
%}
        %% Low Shelf Update & Set Functions

        function update_LS(plugin)
        f0 = plugin.lsFreqSmoother.step();
        Q=plugin.lsQSmoother.step();

        gain = plugin.lsGainSmoother.step();
        mu = db2mag(gain);
        theta = (2 * pi * f0 )/ plugin.fs;
         beta = 4 / (1 + mu);
         delta = beta * tan(theta / 2);
         gamma = (1 - delta) / (1 + delta);

         plugin.filter_LS.a0 = (1 - gamma) / 2;
         plugin.filter_LS.a1 = (1 - gamma) / 2;
         plugin.filter_LS.a2 = 0;
         plugin.filter_LS.b0 = 1;
         plugin.filter_LS.b1 = -gamma;
         plugin.filter_LS.b2 = 0;

%{
        w0=2*pi*f0/plugin.fs;
        alpha=sin(w0)/(2*Q);
        A=sqrt(db2mag(gain));
            
        plugin.filter_LS.a0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha );
        plugin.filter_LS.a1 =  2*A*( (A-1) - (A+1)*cos(w0)                   );
        plugin.filter_LS.a2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha );
        plugin.filter_LS.b0 =          (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha;
        plugin.filter_LS.b1 =   -2*( (A-1) + (A+1)*cos(w0)                   );
        plugin.filter_LS.b2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha;
%}

        end

        function set.LS_FREQ(plugin,val)
        plugin.LS_FREQ = val;
        plugin.lsFreqSmoother.setTargetValue(val);  %#ok  
        update_LS(plugin);
        end

        function set.LS_GAIN(plugin,val)
        plugin.LS_GAIN = val;
        plugin.lsGainSmoother.setTargetValue(val); %#ok
        update_LS(plugin);
        end

        function set.LS_Q(plugin,val)
        plugin.LS_Q = val;    
        plugin.lsQSmoother.setTargetValue(val);
        update_LS(plugin);
        end

       %% High Pass Update & Set functions
        function set.HP_FREQ(plugin,val)
        plugin.HP_FREQ = val;
        plugin.hpFreqSmoother.setTargetValue(val); %#ok  
        update_HP(plugin);
        end

        function set.HP_Q(plugin,val)
        plugin.HP_Q = val;
        plugin.hpQSmoother.setTargetValue(val);
        update_HP(plugin);
        end
        
       function update_HP(plugin)
           f0 = plugin.hpFreqSmoother.step();
           Q = plugin.hpQSmoother.step();

           theta = (2 * pi * f0 )/ plugin.fs;
           d = 1/Q;
           beta = 0.5 * (1 - ((d/2) *(sin(theta)))) / (1 + ((d/2) * (sin(theta))));
           gamma = (0.5 + beta) * cos(theta);

           plugin.filter_HP.a0 = (0.5 + beta + gamma)/2; 
           plugin.filter_HP.a1 = -(0.5 + beta + gamma);
           plugin.filter_HP.a2 = (0.5 + beta + gamma)/2;
           plugin.filter_HP.b0 = 1;
           plugin.filter_HP.b1 = - (2 * gamma); 
           plugin.filter_HP.b2 = 2 * beta;
       end
        %% Low Pass Update & Set functions

        function update_LP(plugin)
           f0 = plugin.lpFreqSmoother.step();
           Q = plugin.lpQSmoother.step();
%{
           w0 = ( 2 * pi * f0)/plugin.fs;

           alpha = sin(w0)/ ( 2 * Q);

           plugin.filter_LP.a0 = (1 - cos(w0))/2;
           plugin.filter_LP.a1 = 1-cos(w0);
           plugin.filter_LP.a2 = (1-cos(w0))/2;
           plugin.filter_LP.b0 = 1 + alpha;
           plugin.filter_LP.b1 = -2 * cos(w0);
           plugin.filter_LP.b2 = 1 - alpha;
        end
%}

           theta = (2 * pi * f0 )/ plugin.fs;
           d = 1/Q;
           beta = 0.5 * (1 - ((d/2) *(sin(theta)))) / (1 + ((d/2) * (sin(theta))));
           gamma = (0.5 + beta) * cos(theta);

           plugin.filter_LP.a0 = (0.5 + beta - gamma)/2; 
           plugin.filter_LP.a1 = (0.5 + beta - gamma);
           plugin.filter_LP.a2 = (0.5 + beta - gamma)/2;
           plugin.filter_LP.b1 = - (2 * gamma); 
           plugin.filter_LP.b2 = 2 * beta;
        end


        function set.LP_FREQ(plugin,val)
        plugin.LP_FREQ = val;
        plugin.lpFreqSmoother.setTargetValue(val); %#ok  
        update_LP(plugin);
        end

        function set.LP_Q(plugin,val)
        plugin.LP_Q = val;
        plugin.lpQSmoother.setTargetValue(val);
        update_LP(plugin);
        end

        %% High Peak Update & Set Functions

        function update_HPF(plugin)
            f0 = plugin.hpfSmoother.step();
            Q = plugin.hpfQSmoother.step();
            gain = plugin.hpfGainSmoother.step();
            %theta = (2 * pi * f0 )/ plugin.fs;
            %mu = db2mag(gain);

            w0 = (2*pi*f0)/plugin.fs;
            alpha = sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));

            plugin.filter_HPF.a0 = 1 + A * alpha;
            plugin.filter_HPF.a1 = -2 * cos(w0);
            plugin.filter_HPF.a2 = 1 - A * alpha;
            plugin.filter_HPF.b0 = 1 + alpha / A;
            plugin.filter_HPF.b1 = -2 * cos(w0);
            plugin.filter_HPF.b2 = 1 - alpha / A;

%{
            zeta = 4 / (1 + mu);
            beta = 0.5 * ((1 - zeta * tan(theta / (2 * Q))) / (1 + zeta * tan(theta / (2 * Q))));
            gamma = (0.5 + beta) * cos(theta);

            plugin.filter_HPF.a0 = 0.5 - beta;
            plugin.filter_HPF.a1 = 0;
            plugin.filter_HPF.a2 = -(0.5 - beta);
            plugin.filter_HPF.b0 = 1;
            plugin.filter_HPF.b1 = -2 * gamma;
            plugin.filter_HPF.b2 = 2 * beta;
%}
        end

        function set.HIGH_PEAK(plugin,val)
        plugin.HIGH_PEAK = val;
        plugin.hpfSmoother.setTargetValue(val);
        update_HPF(plugin);
        end

        function set.HIGHP_GAIN(plugin,val)
        plugin.HIGHP_GAIN = val;
        plugin.hpfGainSmoother.setTargetValue(val);
        update_HPF(plugin);
        end

        function set.HIGHP_Q(plugin,val)
        plugin.HIGHP_Q = val;
        plugin.hpfQSmoother.setTargetValue(val);
        update_HPF(plugin);
        end


        %% Low Peak Update & Set Functions

        function update_LPF(plugin)
            f0 = plugin.lpfSmoother.step();
            Q = plugin.lpfQSmoother.step();
            gain = plugin.lpfGainSmoother.step();

            w0 = 2*pi*f0/plugin.fs;
            alpha = sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));

            plugin.filter_LPF.a0 = 1 + A * alpha;
            plugin.filter_LPF.a1 = -2 * cos(w0);
            plugin.filter_LPF.a2 = 1 - A * alpha;
            plugin.filter_LPF.b0 = 1 + alpha / A;
            plugin.filter_LPF.b1 = -2 * cos(w0);
            plugin.filter_LPF.b2 = 1 - alpha / A;


        end

        function set.LOW_PEAK(plugin,val)
        plugin.LOW_PEAK = val;
        plugin.lpfSmoother.setTargetValue(val);
        update_LPF(plugin);
        end

        function set.LOWP_GAIN(plugin,val)
        plugin.LOWP_GAIN = val;
        plugin.lpfGainSmoother.setTargetValue(val);
        update_LPF(plugin);
        end

        function set.LOWP_Q(plugin,val)
        plugin.LOWP_Q = val;
        plugin.lpfQSmoother.setTargetValue(val);
        update_LPF(plugin);
        end

        %% Mid Peak Update & Set Functions

        function update_MPF(plugin)
            f0 = plugin.mpfSmoother.step();
            Q = plugin.mpfQSmoother.step();
            gain = plugin.mpfGainSmoother.step();

            w0 = 2*pi*f0/plugin.fs;
            alpha = sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));

            plugin.filter_MPF.a0 = 1 + A * alpha;
            plugin.filter_MPF.a1 = -2 * cos(w0);
            plugin.filter_MPF.a2 = 1 - A * alpha;
            plugin.filter_MPF.b0 = 1 + alpha / A;
            plugin.filter_MPF.b1 = -2 * cos(w0);
            plugin.filter_MPF.b2 = 1 - alpha / A;

        end

        function set.MID_PEAK(plugin,val)
        plugin.MID_PEAK = val;
        plugin.mpfSmoother.setTargetValue(val);
        update_MPF(plugin);
        end

        function set.MIDP_GAIN(plugin,val)
        plugin.MIDP_GAIN = val;
        plugin.mpfGainSmoother.setTargetValue(val);
        update_MPF(plugin);
        end

        function set.MIDP_Q(plugin,val)
        plugin.MIDP_Q = val;
        plugin.mpfQSmoother.setTargetValue(val);
        update_MPF(plugin);
        end

        %% Reverb Update & Set Functions

        
        function updateReverb(plugin)
            plugin.reverb.DecayFactor = plugin.DECAY;
            plugin.reverb.WetDryMix = plugin.MIX/100;
        end

        function set.DECAY(plugin,val)
            plugin.DECAY = val;
            updateReverb(plugin);
        end

        function set.MIX(plugin,val)
            plugin.MIX = val;
            updateReverb(plugin);
        end



        %% Gain Update & Set Functions


        function updateGain(plugin)
            plugin.gainLinear = db2mag(plugin.GAIN_dBSmoother.step());
        end

        function set.GAIN_dB(plugin, val)
            plugin.GAIN_dB = val;
            plugin.GAIN_dBSmoother.setTargetValue(val);
            updateGain(plugin);
        end


     end





































end
