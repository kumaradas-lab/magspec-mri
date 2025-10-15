function [t1, T1, data, SeqOut] = get_T1(HW, Seq)
% [t1, T1, data, SeqOut] = get_T1(HW, Seq)
% measures T1 value using sequence_Recovery.m
%
% This function is obsolete and will be removed in a future version of
% OpenMatlab. Please, call "sequence_Recovery" directly.

% Create standard parameters if missing
% if nargin == 1; Seq.Plot=0; end
% if ~isfield(Seq,'tRelax');      Seq.tRelax=[];              end; if isempty(Seq.tRelax);        Seq.tRelax      =   0;                          end
% if ~isfield(Seq,'Flip');        Seq.Flip=[];                end; if isempty(Seq.Flip);          Seq.Flip        =   5;                          end
% if ~isfield(Seq,'nFids');       Seq.nFids=[];               end; if isempty(Seq.nFids);         Seq.nFids       =   50;                         end
% if ~isfield(Seq,'FlipPulse');   Seq.FlipPulse=[];           end; if isempty(Seq.FlipPulse);     Seq.FlipPulse   =   @Pulse_Rect;                end
% if ~isfield(Seq,'InvertPulse'); Seq.InvertPulse=[];         end; if isempty(Seq.InvertPulse);   Seq.InvertPulse =   @Pulse_Rect_Composite180;   end
% if ~isfield(Seq,'Plot');        Seq.Plot=[];                end; if isempty(Seq.Plot);          Seq.Plot        =   0;                          end
% if ~isfield(Seq,'ConsoleOut');  Seq.ConsoleOut=[];          end; if isempty(Seq.ConsoleOut);    Seq.ConsoleOut  =   0;                          end
% if ~isfield(Seq,'tFlip');       Seq.tFlip=[];               end; if isempty(Seq.tFlip);         Seq.tFlip       =   10e-3;                      end
% if ~isfield(Seq,'tFlipStart');  Seq.tFlipStart=[];          end; if isempty(Seq.tFlipStart);    Seq.tFlipStart  =   Seq.tFlip;                  end
% if ~isfield(Seq,'tFlipEnd');    Seq.tFlipEnd=[];            end; if isempty(Seq.tFlipEnd);      Seq.tFlipEnd    =   Seq.tFlip*Seq.nFids;        end
% if ~isfield(Seq,'tFlipLog');    Seq.tFlipLog=[];            end; if isempty(Seq.tFlipLog);      Seq.tFlipLog    =   0;                          end

warning('PD:obsolete_function', ['This function is obsolete and will be removed ' ...
  'in a future version of OpenMatlab. Please, call "sequence_Recovery" directly.']);

% Start measurement by calling the pre-made sequence
[t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq);


%% -----------------------------------------------------------------------------
% (C) Copyright 2012-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
