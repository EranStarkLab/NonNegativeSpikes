% waveform_categorization           categorizes multi-site waveforms
%
% Call: 
%               [ pol, uType, mch, ext, vB ] = waveform_categorization( w, s, THs )
% Gets:
%               w           mean waveform, nsamples x nchannels
%               s           SD of waveform, nsamples x nchannels
%               THs         SD threshold for waveform categorization;
%                               first two values apply for B-spikes, third
%                               value for P-spikes, and forth for N-spikes
%                               {[ 1.25 -1 1.75 -1.75 ]}
%
% Returns:
%               pol         categorization of every waveform (per channel): 
%                               B-spike :=  0
%                               P-spike :=  1
%                               N-spike := -1
%               uType       classification per unit according to the main channel: 
%                               BIP :=      3
%                               Punit :=    2
%                               Other :=    1
%               mch         main channel
%               ext         maximal values of local extrema per channel
%               vB          BPI per channel
% 
% Calls:
%               nothing
%
% Reference:
%               Someck, Levi, Sloin, Spivak, Gattegno, Stark, 2023, 
%               "Positive and biphasic extracellular waveforms correspond to return currents and axonal spikes", 
%               Communications Biology

% 20-aug-23 SS

% last update
% 26-aug-23


function [ pol, uType, mch, ext, vB ] = waveform_categorization( w, s, THs )

% constants
BLS                             = 3;                    % number of samples in the baseline (to subtract from every channel)
USF                             = 4;                    % upsampling factor
TH_B                            = [ -0.6 0.8];          % threshold for BPI definition

% argument handling
nargs                           = nargin;
if nargs < 2 || isempty( w ) || isempty( s )
    error( 'both w and s are required' )
end
[ m, n ]                        = size( w );
if ~isequal( [ m n ], size( s ) )
    error( 'input size mismatch' )
end
BLS                             = round( BLS );
if BLS > m
    BLS                         = m;
elseif BLS < 1
    BLS                         = 1;
end
if USF < 1
    USF                         = 1;
end
USF                             = round( USF );
if nargs < 3 || isempty( THs ) || length( THs( : ) ) ~= 4
    THs                        = [ 1.25 -1 1.75 -1.75 ];
end

% (1) receive mean and SD of waveform in every channel

% (2) subtract the baseline (value in first few samples)
w                               = w - ( ones( m, 1 ) * mean( w( 1 : BLS, : ) ) );

% (3) upsample n-fold 
m2                              = USF * m;
w2                              = zeros( m2, n );
s2                              = zeros( m2, n );
for i                           = 1 : n
    w2( :, i )                  = interp1( 1 : USF : m2, w( :, i ), 1 : m2, 'spline' );
    s2( :, i )                  = interp1( 1 : USF : m2, s( :, i ), 1 : m2, 'spline' );
end
% ignore first/last samples (upsampled) to avoid induction of extrema
nidx                            = [ 1 : USF ( m2 - USF + 1 ) : m2 ];
w2( nidx, : )                   = NaN;
s2( nidx, : )                   = NaN;

% (4) detect all local maxima
[ eidx, evals, etype ]       	= local_find_local_extrema( w2 );

% (5) denote the time, value, and SD of the extremal local minimum/maximum
emat                            = [ eidx evals etype ];
ridx                            = emat( :, 3 ) < 0 & emat( :, 4 ) == 1 | emat( :, 3 ) > 0 & emat( :, 4 ) == -1;
emat( ridx, : )                 = [];                                       % local maximum/minimum must also be positive/negative 
pT                              = NaN( 1, n );                              % sample of positive local extremum per channel                        
pW                              = NaN( 1, n );                              % value of positive local extremum per channel
pS                              = NaN( 1, n );                              % sd for positive local extrema per channel
nT                              = NaN( 1, n );                              % sample of negative local extrema per channel
nW                              = NaN( 1, n );                              % negative local extrema per channel
nS                              = NaN( 1, n );                              % sd for negative local extrema per channel

for i = 1 : n
    ridx                        = emat( :, 2 ) == i;
    imat                        = emat( ridx, : );
    [ ~, maxidx ]               = max( imat( :, 3 ) );
    if imat( maxidx, 4 ) == 1
        pT( i )                 = imat( maxidx, 1 );
        pW( i )                 = w2( imat( maxidx, 1 ), i ); % maxval
        pS( i )                 = s2( imat( maxidx, 1 ), i );
    end
    [ ~, minidx ]               = min( imat( :, 3 ) );
    if imat( minidx, 4 ) == -1
        nT( i )                 = imat( minidx, 1 );
        nW( i )                 = w2( imat( minidx, 1 ), i ); % minval
        nS( i )                 = s2( imat( minidx, 1 ), i );
    end
end
pT                              = pT / USF;
nT                              = nT / USF;

% (6) compute BPI and z-score the extrema
vB                              = ( pW - abs( nW ) ) ./ ( pW + abs( nW ) );
pZ                              = pW ./ pS;
nZ                              = nW ./ nS;
% in case there is no pW/nW, assign -1/1
ex( isnan( pW ) & ~isnan( nW ) ) = -1;
ex( isnan( nW ) & ~isnan( pW ) ) = 1;
vB( ex ~= 0 )                    = ex( ex ~= 0 );

% (7) classification rules
pol                             = NaN ( 1, n );
ext                             = NaN ( 1, n );
for i                           = 1 : n
    if pZ( i ) > THs( 1 ) && nZ( i ) < THs( 2 ) && pT( i ) < nT( i ) && ( vB( i ) > TH_B( 1 ) && vB( i ) < TH_B ( 2 ))
        if ~isnan ( vB( i ) )
            pol( i )            = 0;                                               % B-spike
            ext( i )            = pW(i)-nW(i);
        end
    elseif pZ( i ) > THs( 3 ) && ( pW( i ) > abs( nW( i ) ) || isnan( nW( i ) ) )  % P-spike
        pol( i )                = 1;
        ext( i )                = pW( i );
    elseif nZ( i ) < THs( 4 ) && ( abs( nW( i ) ) > pW( i ) || isnan( pW( i ) ) )  % N-spike
        pol( i )                = -1;
        ext( i )                = nW( i );
    end
end

% (8) unit level categorization
[ ~, mch ] = max( abs( ext ) );                                           	% main channel
upol                            =  pol( mch );
if isempty( upol ) || isnan( upol )
    uType                       = NaN;
elseif upol == 1
    uType                       = 2;
elseif upol == 0
    uType                       = 3;
elseif upol == -1 
    uType                       = 1;    
end

return % waveform_categorization.m

%------------------------------------------------------------------------
% [ idx, vals, etype ] = local_find_local_extrema( x )
% detect all local extrema in a matrix
%------------------------------------------------------------------------
function [ idx, vals, etype ] = local_find_local_extrema( x )

if isempty( x )
    return
end
sx                              = size( x );
if length( sx ) == 2 && sum( sx == 1 ) >= 1 && sx( 2 ) > 1
    x                           = x';
end
if length( sx ) > 2
    error( 'input size mismatch: x must be a vector or a matrix' )
end
m                               = sx( 1 );
n                               = sx( 2 );

% compute second 'derivative'
d2                              = diff( sign( diff( x ) ) );

% identify extrema
[ rowMin, colMin ]              = find( d2 > 1 );
[ rowMax, colMax ]              = find( d2 < -1 );
etype                           = [ -1 * ones( length( rowMin ), 1 ); ones( length( rowMax ), 1 ) ];
mat                             = [ [ [ rowMin colMin ]; [ rowMax colMax ] ] etype ];
smat                            = sortrows( mat, [ 2 1 ] );
row                             = smat( :, 1 );
col                             = smat( :, 2 );
etype                           = smat( :, 3 );
row                             = row + 1;

% organize output
if n == 1
    idx                         = row;
else
    idx                         = [ row col ];
end
vals                            = x( row + ( col - 1 ) * m );

return % local_find_local_extrema

%------------------------------------------------------------------------
% EOF
%------------------------------------------------------------------------

s_file                          = 's_exampledata_mDS2_07.mat';              % example dataset (137 units from a single neocortex session)
%s_file                          = 's_full_dataset.mat';                     % full dataset (9160 units from 197 sessions)
load( s_file, '-mat', 's' )

n                               = size( s.mean, 1 );
nsites_all                      = NaN( n, 1 );
for i                           = 1 : n
    nsites_all( i )             = size( s.mean{ i }, 1 );
end
nsites                          = max( nsites_all );
z0                              = NaN( n, nsites );
pol                             = z0; 
uType                           = NaN( n, 1 );
mch                             = NaN( n, 1 );
ext                             = z0;
vB                              = z0;

% go over units and compute stats
for i                           = 1 : n 
    w                           = s.mean{ i };
    sd                          = s.sd{ i };
    [pol0, uType0, mch0 ...
        , ext0, vB0]            = waveform_categorization( w', sd' );
    % update the structure
    cidx                        = 1 : length( pol0 );
    pol( i, cidx )              = pol0;
    uType( i )                  = uType0;
    mch( i )                  	= mch0;
    ext( i, cidx )           	= ext0;
    vB( i, cidx )           	= vB0;
end

% expand the structure with the new fields
s.pol                           = pol;
s.uType                         = uType;
s.mch                           = mch;
s.ext                           = ext;
s.vB                            = vB;
