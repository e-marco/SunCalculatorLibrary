{
Copyright (c) 2015 Institute for Solar Energy Research Hamelin
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
 4. Redistributions of any form whatsoever must retain the following
    acknowledgment: 'This product includes the SunCalculator-Library developed 
	by Marco Ernst and Hendrik Holst at the 
	Institute for Solar Energy Research Hamelin (www.isfh.de)'
	
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
}
library binning;

{ Important note about DLL memory management: ShareMem must be the
  first unit in your library's USES clause AND your project's (select
  Project-View Source) USES clause if your DLL exports any procedures or
  functions that pass strings as parameters or function results. This
  applies to all strings passed to and from your DLL--even those that
  are nested in records and classes. ShareMem is the interface unit to
  the BORLNDMM.DLL shared memory manager, which must be deployed along
  with your DLL. To avoid using BORLNDMM.DLL, pass string information
  using PChar or ShortString parameters. }

uses
  System.SimpleShareMem,
  System.SysUtils,
  System.Classes,
  Vcl.Grids,
  System.Math,
  radiationTypes in 'radiationTypes.pas';

type TArrayOfDouble = Array of double;

var
  BinMatrix : array [1..752,1..129600] of double;
  PostProcessingMatrix : array [1..752,1..129600] of double;

  isBinMatrixInitialized : boolean;
  isPostProcessingMatrixInitialized : boolean;

const
  wantPhotonFlux = 0; // 0 = no, 1 = yes

{$R *.res}

function resetInitialization:boolean; STDCALL;
begin
  isPostProcessingMatrixInitialized := false;
  isBinMatrixInitialized := False;
end;


// Call before first input of binning data
function InitializeBinning(BinWidth:integer;coordinateSystem:integer):boolean; STDCALL;
var i,j,x,y : integer;
    x_max, y_max : integer;
    max_i : integer;
begin

  try
    result := true;
    isBinMatrixInitialized := True;

    max_i := ceil((360 * 360) / (BinWidth * BinWidth));

    for I := 1 to max_i do
      for j := 3 to 752 do
        BinMatrix[j,i] := 0;  // clear z Axis

    if (coordinateSystem=0) or (coordinateSystem=3) then // horizontal and panel coordinates
    begin
      x_max:=360; // AZI = 0° ... 360°
      y_max:=180; // ALT = -90° ... 90°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        y := 1;
          while y<=y_max do
          begin
            BinMatrix[1,i] := (x-1) + BinWidth/2;
            BinMatrix[2,i] := (y-1) + BinWidth/2 - 90;
            y := y + BinWidth;
            i := i + 1;
          end;
        x := x + BinWidth;
      end;
    end;

    if coordinateSystem=1 then
    begin
      x_max:=180; // YAN = -90° ... 90°
      y_max:=180; // ZAN = -90° ... 90°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        y := 1;
          while y<=y_max do
          begin
            BinMatrix[1,i] := (x-1) + BinWidth/2 - 90;
            BinMatrix[2,i] := (y-1) + BinWidth/2 - 90;
            y := y + BinWidth;
            i := i + 1;
          end;
        x := x + BinWidth;
      end;
    end;

    if coordinateSystem=2 then
    begin
      x_max:=180; // TTA = 0° ... 180°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        BinMatrix[1,i] := (x-1) + BinWidth/2;
        x := x + BinWidth;
        i:=i+1;
      end;
    end;

  except
    result := false;
  end;
end;

// Call before first input of post processing data
function InitializePostProcessing(BinWidth:integer;coordinateSystem:integer):boolean; STDCALL;
var i,j,x,y : integer;
    x_max, y_max, max_i : integer;
begin
  try
    result := true;
    isPostProcessingMatrixInitialized := True;

    max_i := ceil((360 * 360) / (BinWidth * BinWidth));

    for I := 1 to max_i do
      for j := 3 to 752 do
        PostProcessingMatrix[j,i] := 0;  // clear z Axis

    if (coordinateSystem=0) or (coordinateSystem=3) then
    begin
      x_max:=360; // AZI = 0° ... 360°
      y_max:=180; // ALT = -90° ... 90°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        y := 1;
          while y<=y_max do
          begin
            PostProcessingMatrix[1,i] := (x-1) + BinWidth/2;
            PostProcessingMatrix[2,i] := (y-1) + BinWidth/2 - 90;
            y := y + BinWidth;
            i := i + 1;
          end;
        x := x + BinWidth;
      end;

    end;

    if coordinateSystem=1 then
    begin
      x_max:=180; // YAN = -90° ... 90°
      y_max:=180; // ZAN = -90° ... 90°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        y := 1;
          while y<=y_max do
          begin
            PostProcessingMatrix[1,i] := (x-1) + BinWidth/2 - 90;
            PostProcessingMatrix[2,i] := (y-1) + BinWidth/2 - 90;
            y := y + BinWidth;
            i := i + 1;
          end;
        x := x + BinWidth;
      end;
    end;

    if coordinateSystem=2 then
    begin
      x_max:=180; // TTA = 0° ... 180°

      x:=1;
      i:=1;
      while x<=x_max do
      begin
        PostProcessingMatrix[1,i] := (x-1) + BinWidth/2;
        x := x + BinWidth;
        i:=i+1;
      end;
    end;

  except
    result := false;
  end;

end;

// Bin photon flux data either by ALT and AZI or by YAN and ZAN.
// coordinateSystem=0 -> the function uses ALT & AZI, values of YAN & ZAN are ignored
// ALT_AZI=false -> the function uses YAN & ZAN, values of ALT & AZI are ignored
// z_value is the z-data, could be any floating point value
// Direct is true for direct irradiance and false for diffuse irradiance
function BinValues(AZI_value,ALT_value,YAN_value,ZAN_value,TTA_value,irradiance,photonflux:double;BinWidth:integer;
  coordinateSystem:integer; IsDirect:boolean; wantPhotonData:boolean; photonData:TArrayOfDouble):boolean; STDCALL;
var i,counter,numberOfWavelengths : integer;
    x_max, y_max, x_bins, y_bins : integer;
    x_binwidth, y_binwidth : integer;
    x_val, y_val : double;
    x_search, y_search:double;
    x_row, y_row, findrow : integer;
label
  RowLocatedAA,RowLocatedYZ,RowLocatedTTA;
begin
  try
    result := true;
    if wantPhotonData then
      numberOfWavelengths := Length(photonData)
    else
      numberOfWavelengths := 0;

    if (coordinateSystem=0) or (coordinateSystem=3) then
    begin
      y_max:=180; // ALT = -90° ... 90°

    // Suche Reihe
      AZI_value := min(AZI_value, 359.999);

      findRow:=Floor(AZI_value/binWidth)*Round(y_max/BinWidth) + Floor((ALT_value+89.999)/BinWidth) + 1;


      RowLocatedAA:

      if findrow=-1 then
        result := false else
          if IsDirect then
          begin     // direct or diffuse radiation?
            BinMatrix[3, findrow] := BinMatrix[3, findrow] + irradiance;
            BinMatrix[5, findrow] := BinMatrix[5, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + counter,findRow]:=BinMatrix[7 + counter,findRow] + photonData[counter];
              end;
            end;
          end
          else
          begin
            BinMatrix[4, findrow] := BinMatrix[4, findrow] + irradiance;
            BinMatrix[6, findrow] := BinMatrix[6, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + numberOfWavelengths + counter,findRow]:=BinMatrix[7 + numberOfWavelengths + counter,findRow] + photonData[counter];
              end;
            end;
          end;
    end;

    if coordinateSystem=1 then
    begin
      x_max:=180; // YAN = -90° ... 90°
      y_max:=180; // ZAN = -90° ... 90°

      x_binwidth:=BinWidth;
      y_binwidth:=BinWidth;

      x_bins:=ceil(x_max / x_binwidth);
      y_bins:=ceil(y_max / y_binwidth);

      // Suche Reihe

      x_val := Floor(YAN_value + 90);
      x_search := (Floor(x_val) div BinWidth)*BinWidth - 90 + BinWidth/2;

      y_val := Floor(ZAN_value + 90);
      y_search := (Floor(y_val) div BinWidth)*BinWidth - 90 + BinWidth/2;

      findrow:=-1;

      for x_row := 1 to x_bins do
        for y_row := 1 to y_bins do
        begin
          i := y_row + (x_row-1)*y_bins;
          if (x_search=BinMatrix[1,i]) and (y_search=BinMatrix[2,i])
            then begin
              findrow := i;
              Goto RowLocatedYZ;
            end;
        end;

      RowLocatedYZ:

      if findrow=-1 then
        result := false else
          if IsDirect then
          begin     // direct or diffuse radiation?
            BinMatrix[3, findrow] := BinMatrix[3, findrow] + irradiance;
            BinMatrix[5, findrow] := BinMatrix[5, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + counter,findRow]:=BinMatrix[7+ counter,findRow] + photonData[counter];
              end;
            end;
          end
          else
          begin
            BinMatrix[4, findrow] := BinMatrix[4, findrow] + irradiance;
            BinMatrix[6, findrow] := BinMatrix[6, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + numberOfWavelengths + counter,findRow]:=BinMatrix[7 + numberOfWavelengths + counter,findRow] + photonData[counter];
              end;
            end;
          end;
    end;

    if coordinateSystem=2 then
    begin
      x_max:=180; // TTA = 0° ... 180°

      x_binwidth:=BinWidth;

      x_bins:=ceil(x_max / x_binwidth);

      // Suche Reihe

      x_search := (Floor(TTA_value) div BinWidth)*BinWidth + BinWidth/2;

      findrow:=-1;

      for x_row := 1 to x_bins do
        if (x_search=BinMatrix[1,x_row]) then
        begin
          findrow := x_row;
          Goto RowLocatedTTA;
        end;

      RowLocatedTTA:

      if findrow=-1 then
        result := false else
          if IsDirect then
          begin     // direct or diffuse radiation?
            BinMatrix[3, findrow] := BinMatrix[3, findrow] + irradiance;
            BinMatrix[5, findrow] := BinMatrix[5, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + counter,findRow]:=BinMatrix[7+ counter,findRow] + photonData[counter];
              end;
            end;
          end
          else
          begin
            BinMatrix[4, findrow] := BinMatrix[4, findrow] + irradiance;
            BinMatrix[6, findrow] := BinMatrix[6, findrow] + photonflux;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths - 1 do
              begin
                BinMatrix[7 + numberOfWavelengths + counter,findRow]:=BinMatrix[7 + numberOfWavelengths + counter,findRow] + photonData[counter];
              end;
            end;
          end;
    end;
  except
    result := false;
  end;
end;


// Call to get the binning data
// coordinateSystem=0 -> the function uses ALT & AZI, values of YAN & ZAN are ignored
// ALT_AZI=false -> the function uses YAN & ZAN, values of ALT & AZI are ignored
// PhotonFlux is the z-data, could be any floating point value
// direct, diffuse : true if checked!
function BinningOutput(BinWidth : integer; coordinateSystem:integer;
  DirectIrradiance, DiffuseIrradiance:boolean; irradiance_title, flux_title:string;
  wantPhotonDataDirect, wantPhotonDataDiffuse:boolean; numberOfWavelengths:integer; startWavelength, endWavelength, intervalWavelength : integer):TStringGrid; STDCALL;
var j,counter,photonCount : integer;
    x_max, y_max, x_bins, y_bins : integer;
    x_binwidth, y_binwidth : integer;
begin
  result := TStringGrid.Create(nil);
  if(wantPhotonDataDirect) and (wantPhotonDataDiffuse) then
    photonCount := 2 * numberOfWavelengths
  else if (wantPhotonDataDirect) or (wantPhotonDataDiffuse) then
    photonCount := numberOfWavelengths
  else
    photonCount := 0;


  if (coordinateSystem=0) or (coordinateSystem=3) then
  begin
    x_max:=360; // AZI = 0° ... 360°
    y_max:=180; // ALT = -90° ... 90°

    x_binwidth:=BinWidth;
    y_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
    y_bins:=ceil(y_max / y_binwidth);
  end;

  if coordinateSystem=1 then
  begin
    x_max:=180; // YAN = -90° ... 90°
    y_max:=180; // ZAN = -90° ... 90°

    x_binwidth:=BinWidth;
    y_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
    y_bins:=ceil(y_max / y_binwidth);
  end;

  if coordinateSystem=2 then
  begin
    x_max:=180; // TTA = 0° ... 180°

    x_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
  end;


  with result do
  begin
    if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
    begin
      RowCount:=x_bins*y_bins+1;      // +1 for title
      if DirectIrradiance and DiffuseIrradiance then
        ColCount := 4 + 2*wantPhotonFlux + photonCount
      else
        ColCount := 3 + wantPhotonFlux + photonCount;
    end;

    if (coordinateSystem=2) then
    begin
      RowCount:=x_bins+1;      // +1 for title
      if DirectIrradiance and DiffuseIrradiance then
        ColCount := 3 + 2*wantPhotonFlux + photonCount
      else
        ColCount := 2 + wantPhotonFlux + photonCount;
    end;

    // title

    if (coordinateSystem=0) or (coordinateSystem=3) then
    begin
      if (coordinateSystem=0) then
      begin
        Cells[0,0]:='Azimuth_(AZI)';
        Cells[1,0]:='Altitude_(ALT)';
      end else begin
        Cells[0,0]:='PanelAzimuth_(AZI)';
        Cells[1,0]:='PanelAltitude_(ALT)';
      end;

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[2,0]:='Direct '+irradiance_title;
        Cells[3,0]:='Diffuse '+irradiance_title;
        if wantPhotonFlux=1 then
        begin
          Cells[4,0]:='Direct '+flux_title;
          Cells[5,0]:='Diffuse '+flux_title;
        end;
        if wantPhotonDataDirect then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonDataDiffuse then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j+photonCount - numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;

      end else
        if DirectIrradiance then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          if wantPhotonDataDirect then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

          if wantPhotonflux=1 then
            Cells[3,0]:='Direct '+flux_title;
        end
        else begin
          Cells[2,0]:='Diffuse '+irradiance_title;
          if wantPhotonDataDiffuse then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          if wantPhotonflux=1 then
            Cells[3,0]:='Diffuse '+flux_title;
        end;


    end;

    if coordinateSystem=1 then
    begin
      Cells[0,0]:='Y-Angle_(YAN)';
      Cells[1,0]:='Z-Angle_(ZAN)';

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[2,0]:='Direct '+irradiance_title;
        Cells[3,0]:='Diffuse '+irradiance_title;
        if wantPhotonDataDirect then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;

        if wantPhotonDataDiffuse then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j+photoncount - numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonflux=1 then
        begin
          Cells[4,0]:='Direct '+flux_title;
          Cells[5,0]:='Diffuse '+flux_title;
        end;
      end else
        if DirectIrradiance then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[3,0]:='Direct '+flux_title;
          if wantPhotonDataDirect then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

        end
        else begin
          Cells[2,0]:='Diffuse '+irradiance_title;
          if wantPhotonDataDiffuse then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          if wantPhotonflux=1 then
            Cells[3,0]:='Diffuse '+flux_title;
        end;
    end;

    if coordinateSystem=2 then
    begin
      Cells[0,0]:='Inclination-Angle_(TTA)';

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[1,0]:='Direct '+irradiance_title;
        Cells[2,0]:='Diffuse '+irradiance_title;
        if wantPhotonDataDirect then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[3+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonDataDiffuse then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[3+2*wantPhotonflux+j+photoncount - numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonflux=1 then
        begin
          Cells[3,0]:='Direct '+flux_title;
          Cells[4,0]:='Diffuse '+flux_title;
        end;
      end else
        if DirectIrradiance then
        begin
          Cells[1,0]:='Direct '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[2,0]:='Direct '+flux_title;
          if wantPhotonDataDirect then
            for j := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

        end
        else begin
          Cells[1,0]:='Diffuse '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[2,0]:='Diffuse '+flux_title;
          if wantPhotonDataDiffuse then
            for j := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
    end;

    // data

    if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
      for j := 1 to RowCount do
      begin
        Cells[0,j]:=FloatToStr(BinMatrix[1,j]); // AZI bzw. YAN
        Cells[1,j]:=FloatToStr(BinMatrix[2,j]); // ALT bzw. ZAN

        if DirectIrradiance and DiffuseIrradiance then
        begin
          Cells[2,j]:=FloatToStrF(BinMatrix[3,j], ffGeneral, 8,1); // Direct irradiance
          Cells[3,j]:=FloatToStrF(BinMatrix[4,j], ffGeneral, 8,1); // Diffuse irradiance
          if wantPhotonflux=1 then
          begin
            Cells[4,j]:=FloatToStrF(BinMatrix[5,j], ffFixed, 13,0); // Direct photon flux
            Cells[5,j]:=FloatToStrF(BinMatrix[6,j], ffFixed, 13,0); // Diffuse photon flux
          end;
          if wantPhotonDataDirect then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+counter,j]:=FloatToStrF(BinMatrix[7+counter,j], ffFixed, 13,0);
          end;
          if wantPhotonDataDiffuse then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+counter+photoncount - numberOfWavelengths,j]:=FloatToStrF(BinMatrix[7+numberOfWavelengths+counter,j], ffFixed, 13,0);
          end;
          
        end else
        if DirectIrradiance then
        begin
          Cells[2,j]:=FloatToStrF(BinMatrix[3,j], ffGeneral, 8,1); // Direct irradiance
          if wantPhotonflux=1 then
            Cells[3,j]:=FloatToStrF(BinMatrix[5,j], ffFixed, 15,0); // Direct photon flux
          if wantPhotonDataDirect then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+counter,j]:=FloatToStrF(BinMatrix[7+counter,j], ffFixed, 13,0);
        end
        else begin
          Cells[2,j]:=FloatToStrF(BinMatrix[4,j], ffGeneral, 8,1); // Diffuse irradiance
          if wantPhotonflux=1 then
            Cells[3,j]:=FloatToStrF(BinMatrix[6,j], ffFixed, 15,0); // Diffuse photon flux
          if wantPhotonDataDiffuse then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+counter+photoncount - numberOfWavelengths,j]:=FloatToStrF(BinMatrix[7+numberOfWavelengths+counter,j], ffFixed, 13,0);
        end;
      end;

    if (coordinateSystem=2) then
      for j := 1 to RowCount do
      begin
        Cells[0,j]:=FloatToStr(BinMatrix[1,j]); // Matrix export

        if DirectIrradiance and DiffuseIrradiance then
        begin
          Cells[1,j]:=FloatToStrF(BinMatrix[3,j], ffGeneral, 8,1); // Direct irradiance
          Cells[2,j]:=FloatToStrF(BinMatrix[4,j], ffGeneral, 8,1); // Diffuse irradiance
          if wantPhotonflux=1 then
          begin
            Cells[3,j]:=FloatToStrF(BinMatrix[5,j], ffFixed, 15,0); // Direct photon flux
            Cells[4,j]:=FloatToStrF(BinMatrix[6,j], ffFixed, 15,0); // Diffuse photon flux
          end;
          if wantPhotonDataDirect then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+counter,j]:=FloatToStrF(BinMatrix[7+counter,j], ffFixed, 13,0);
          end;
          if wantPhotonDataDiffuse then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+counter+photoncount - numberOfWavelengths,j]:=FloatToStrF(BinMatrix[7+numberOfWavelengths+counter,j], ffFixed, 13,0);
          end;
        end else
        if DirectIrradiance then
        begin
          Cells[1,j]:=FloatToStrF(BinMatrix[3,j], ffGeneral, 8,1); // Direct irradiance
          if wantPhotonflux=1 then
            Cells[2,j]:=FloatToStrF(BinMatrix[5,j], ffFixed, 15,0); // Direct photon flux
          if wantPhotonDataDirect then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+counter,j]:=FloatToStrF(BinMatrix[7+counter,j], ffFixed, 13,0);
        end
        else begin
          Cells[1,j]:=FloatToStrF(BinMatrix[4,j], ffGeneral, 8,1); // Diffuse irradiance
          if wantPhotonflux=1 then
            Cells[2,j]:=FloatToStrF(BinMatrix[6,j], ffFixed, 15,0); // Diffuse photon flux
          if wantPhotonDataDiffuse then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+counter+photoncount-numberOfWavelengths,j]:=FloatToStrF(BinMatrix[7+numberOfWavelengths+counter,j], ffFixed, 13,0);
        end;
      end;
  end;  // with OutputGrid do begin
end;


// Call to get the binning data
// coordinateSystem=0 -> the function uses ALT & AZI, values of YAN & ZAN are ignored
// ALT_AZI=false -> the function uses YAN & ZAN, values of ALT & AZI are ignored
// PhotonFlux is the z-data, could be any floating point value
// direct, diffuse : true if checked!
function BinningOutputVar(BinWidth : integer; coordinateSystem:integer;
  DirectIrradiance, DiffuseIrradiance:boolean; irradiance_title, flux_title:string;
  wantPhotonData:boolean; numberOfWavelengths:integer; startWavelength, endWavelength, intervalWavelength : integer; var OutputGrid : TStringGrid):boolean; STDCALL;
var j,counter,photonCount : integer;
    x_max, y_max, x_bins, y_bins : integer;
    x_binwidth, y_binwidth : integer;
begin
  try
      result := true;

    if(wantPhotonData) then
      photonCount := numberOfWavelengths
    else
      photonCount := 0;


    if (coordinateSystem=0) or (coordinateSystem=3) then
    begin
      x_max:=360; // AZI = 0° ... 360°
      y_max:=180; // ALT = -90° ... 90°

      x_binwidth:=BinWidth;
      y_binwidth:=BinWidth;

      x_bins:=ceil(x_max / x_binwidth);
      y_bins:=ceil(y_max / y_binwidth);
    end;

    if coordinateSystem=1 then
    begin
      x_max:=180; // YAN = -90° ... 90°
      y_max:=180; // ZAN = -90° ... 90°

      x_binwidth:=BinWidth;
      y_binwidth:=BinWidth;

      x_bins:=ceil(x_max / x_binwidth);
      y_bins:=ceil(y_max / y_binwidth);
    end;

    if coordinateSystem=2 then
    begin
      x_max:=180; // TTA = 0° ... 180°

      x_binwidth:=BinWidth;

      x_bins:=ceil(x_max / x_binwidth);
    end;


    with OutputGrid do
    begin
      if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
      begin
        RowCount:=x_bins*y_bins+1;      // +1 for title
        if DirectIrradiance and DiffuseIrradiance then
          ColCount := 4 + 2*wantPhotonFlux + 2*photonCount
        else
          ColCount := 3 + wantPhotonFlux + 2*photonCount;
      end;

      if (coordinateSystem=2) then
      begin
        RowCount:=x_bins+1;      // +1 for title
        if DirectIrradiance and DiffuseIrradiance then
          ColCount := 3 + 2*wantPhotonFlux + 2*photonCount
        else
          ColCount := 2 + wantPhotonFlux + 2*photonCount;
      end;

      // title

      if (coordinateSystem=0) or (coordinateSystem=3) then
      begin
        if (coordinateSystem=0) then begin
          Cells[0,0]:='Azimuth_(AZI)';
          Cells[1,0]:='Altitude_(ALT)';
        end else begin
          Cells[0,0]:='PanelAzimuth_(AZI)';
          Cells[1,0]:='PanelAltitude_(ALT)';
        end;

        if (DirectIrradiance) and (DiffuseIrradiance) then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          Cells[3,0]:='Diffuse '+irradiance_title;
          if wantPhotonFlux=1 then
          begin
            Cells[4,0]:='Direct '+flux_title;
            Cells[5,0]:='Diffuse '+flux_title;
          end;
          if wantPhotonData then
          begin
            for j := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
            for j := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          end;

        end else
          if DirectIrradiance then
          begin
            Cells[2,0]:='Direct '+irradiance_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

            if wantPhotonflux=1 then
              Cells[3,0]:='Direct '+flux_title;
          end
          else begin
            Cells[2,0]:='Diffuse '+irradiance_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
            if wantPhotonflux=1 then
              Cells[3,0]:='Diffuse '+flux_title;
          end;


      end;

      if coordinateSystem=1 then
      begin
        Cells[0,0]:='Y-Angle_(YAN)';
        Cells[1,0]:='Z-Angle_(ZAN)';

        if (DirectIrradiance) and (DiffuseIrradiance) then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          Cells[3,0]:='Diffuse '+irradiance_title;
          if wantPhotonData then
          begin
            for j := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
            for j := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          end;
          if wantPhotonflux=1 then
          begin
            Cells[4,0]:='Direct '+flux_title;
            Cells[5,0]:='Diffuse '+flux_title;
          end;
        end else
          if DirectIrradiance then
          begin
            Cells[2,0]:='Direct '+irradiance_title;
            if wantPhotonflux=1 then
              Cells[3,0]:='Direct '+flux_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

          end
          else begin
            Cells[2,0]:='Diffuse '+irradiance_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
            if wantPhotonflux=1 then
              Cells[3,0]:='Diffuse '+flux_title;
          end;
      end;

      if coordinateSystem=2 then
      begin
        Cells[0,0]:='Inclination-Angle_(TTA)';

        if (DirectIrradiance) and (DiffuseIrradiance) then
        begin
          Cells[1,0]:='Direct '+irradiance_title;
          Cells[2,0]:='Diffuse '+irradiance_title;
          if wantPhotonData then
          begin
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          end;
          if wantPhotonflux=1 then
          begin
            Cells[3,0]:='Direct '+flux_title;
            Cells[4,0]:='Diffuse '+flux_title;
          end;
        end else
          if DirectIrradiance then
          begin
            Cells[1,0]:='Direct '+irradiance_title;
            if wantPhotonflux=1 then
              Cells[2,0]:='Direct '+flux_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[2+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

          end
          else begin
            Cells[1,0]:='Diffuse '+irradiance_title;
            if wantPhotonflux=1 then
              Cells[2,0]:='Diffuse '+flux_title;
            if wantPhotonData then
              for j := 0 to numberOfWavelengths-1 do
                Cells[2+wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          end;
      end;

      // data

      if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
        for j := 1 to RowCount do
        begin
          Cells[0,j]:=FloatToStr(BinMatrix[1,j]); // AZI bzw. YAN
          Cells[1,j]:=FloatToStr(BinMatrix[2,j]); // ALT bzw. ZAN

          if DirectIrradiance and DiffuseIrradiance then
          begin
            Cells[2,j]:=FloatToStr(BinMatrix[3,j]); // Direct irradiance
            Cells[3,j]:=FloatToStr(BinMatrix[4,j]); // Diffuse irradiance
            if wantPhotonflux=1 then
            begin
              Cells[4,j]:=FloatToStr(BinMatrix[5,j]); // Direct photon flux
              Cells[5,j]:=FloatToStr(BinMatrix[6,j]); // Diffuse photon flux
            end;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths-1 do
                Cells[4+2*wantPhotonflux+counter,j]:=FloatToStr(BinMatrix[7+counter,j]);
              for counter := 0 to numberOfWavelengths-1 do
                Cells[4+2*wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(BinMatrix[7+numberOfWavelengths+counter,j]);
            end;

          end else
          if DirectIrradiance then
          begin
            Cells[2,j]:=FloatToStr(BinMatrix[3,j]); // Direct irradiance
            if wantPhotonflux=1 then
              Cells[3,j]:=FloatToStr(BinMatrix[5,j]); // Direct photon flux
            if wantPhotonData then
              for counter := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+counter,j]:=FloatToStr(BinMatrix[7+counter,j]);
          end
          else begin
            Cells[2,j]:=FloatToStr(BinMatrix[4,j]); // Diffuse irradiance
            if wantPhotonflux=1 then
              Cells[3,j]:=FloatToStr(BinMatrix[6,j]); // Diffuse photon flux
            if wantPhotonData then
              for counter := 0 to numberOfWavelengths-1 do
                Cells[3+wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(BinMatrix[7+numberOfWavelengths+counter,j]);
          end;
        end;

      if (coordinateSystem=2) then
        for j := 1 to RowCount do
        begin
          Cells[0,j]:=FloatToStr(BinMatrix[1,j]); // Matrix export

          if DirectIrradiance and DiffuseIrradiance then
          begin
            Cells[1,j]:=FloatToStr(BinMatrix[3,j]); // Direct irradiance
            Cells[2,j]:=FloatToStr(BinMatrix[4,j]); // Diffuse irradiance
            if wantPhotonflux=1 then
            begin
              Cells[3,j]:=FloatToStr(BinMatrix[5,j]); // Direct photon flux
              Cells[4,j]:=FloatToStr(BinMatrix[6,j]); // Diffuse photon flux
            end;
            if wantPhotonData then
            begin
              for counter := 0 to numberOfWavelengths-1 do
                Cells[3+2*wantPhotonflux+counter,j]:=FloatToStr(BinMatrix[7+counter,j]);
              for counter := 0 to numberOfWavelengths-1 do
                Cells[3+2*wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(BinMatrix[7+numberOfWavelengths+counter,j]);
            end;
          end else
          if DirectIrradiance then
          begin
            Cells[1,j]:=FloatToStr(BinMatrix[3,j]); // Direct irradiance
            if wantPhotonflux=1 then
              Cells[2,j]:=FloatToStr(BinMatrix[5,j]); // Direct photon flux
            if wantPhotonData then
              for counter := 0 to numberOfWavelengths-1 do
                Cells[2+wantPhotonflux+counter,j]:=FloatToStr(BinMatrix[7+counter,j]);
          end
          else begin
            Cells[1,j]:=FloatToStr(BinMatrix[4,j]); // Diffuse irradiance
            if wantPhotonflux=1 then
              Cells[2,j]:=FloatToStr(BinMatrix[6,j]); // Diffuse photon flux
            if wantPhotonData then
              for counter := 0 to numberOfWavelengths-1 do
                Cells[2+wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(BinMatrix[7+numberOfWavelengths+counter,j]);
          end;
        end;
    end;  // with OutputGrid do begin
  except
    result := false;
  end;
end;

// Call to get the post processing data
// coordinateSystem=0 -> the function uses ALT & AZI, values of YAN & ZAN are ignored
// ALT_AZI=false -> the function uses YAN & ZAN, values of ALT & AZI are ignored
// PhotonFlux is the z-data, could be any floating point value
// direct, diffuse : true if checked!
function PostProcessingOutput(BinWidth : integer; coordinateSystem:integer;
  DirectIrradiance, DiffuseIrradiance:boolean; irradiance_title, flux_title:string;
  wantPhotonData:boolean; numberOfWavelengths:integer; startWavelength, endWavelength, intervalWavelength : integer):TStringGrid; STDCALL;
var j,counter,photonCount : integer;
    x_max, y_max, x_bins, y_bins : integer;
    x_binwidth, y_binwidth : integer;
begin
 result := TStringGrid.Create(nil);
  if(wantPhotonData) then
    photonCount := numberOfWavelengths
  else
    photonCount := 0;


  if (coordinateSystem=0) or (coordinateSystem=3) then
  begin
    x_max:=360; // AZI = 0° ... 360°
    y_max:=180; // ALT = -90° ... 90°

    x_binwidth:=BinWidth;
    y_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
    y_bins:=ceil(y_max / y_binwidth);
  end;

  if coordinateSystem=1 then
  begin
    x_max:=180; // YAN = -90° ... 90°
    y_max:=180; // ZAN = -90° ... 90°

    x_binwidth:=BinWidth;
    y_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
    y_bins:=ceil(y_max / y_binwidth);
  end;

  if coordinateSystem=2 then
  begin
    x_max:=180; // TTA = 0° ... 180°

    x_binwidth:=BinWidth;

    x_bins:=ceil(x_max / x_binwidth);
  end;


  with result do
  begin
    if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
    begin
      RowCount:=x_bins*y_bins+1;      // +1 for title
      if DirectIrradiance and DiffuseIrradiance then
        ColCount := 4 + 2*wantPhotonFlux + 2*photonCount
      else
        ColCount := 3 + wantPhotonFlux + 2*photonCount;
    end;

    if (coordinateSystem=2) then
    begin
      RowCount:=x_bins+1;      // +1 for title
      if DirectIrradiance and DiffuseIrradiance then
        ColCount := 1 + 2*wantPhotonFlux + 2*photonCount
      else
        ColCount := 1 + wantPhotonFlux + 2*photonCount;
    end;

    // title

    if (coordinateSystem=0) or (coordinateSystem=3) then
    begin
      if (coordinateSystem=0) then begin
        Cells[0,0]:='Azimuth_(AZI)';
        Cells[1,0]:='Altitude_(ALT)';
      end else begin
        Cells[0,0]:='PanelAzimuth_(AZI)';
        Cells[1,0]:='PanelAltitude_(ALT)';
      end;

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[2,0]:='Direct '+irradiance_title;
        Cells[3,0]:='Diffuse '+irradiance_title;
        if wantPhotonFlux=1 then
        begin
          Cells[4,0]:='Direct '+flux_title;
          Cells[5,0]:='Diffuse '+flux_title;
        end;
        if wantPhotonData then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;

      end else
        if DirectIrradiance then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

          if wantPhotonflux=1 then
            Cells[3,0]:='Direct '+flux_title;
        end
        else begin
          Cells[2,0]:='Diffuse '+irradiance_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          if wantPhotonflux=1 then
            Cells[3,0]:='Diffuse '+flux_title;
        end;


    end;

    if coordinateSystem=1 then
    begin
      Cells[0,0]:='Y-Angle_(YAN)';
      Cells[1,0]:='Z-Angle_(ZAN)';

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[2,0]:='Direct '+irradiance_title;
        Cells[3,0]:='Diffuse '+irradiance_title;
        if wantPhotonData then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          for j := 0 to numberOfWavelengths-1 do
            Cells[4+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonflux=1 then
        begin
          Cells[4,0]:='Direct '+flux_title;
          Cells[5,0]:='Diffuse '+flux_title;
        end;
      end else
        if DirectIrradiance then
        begin
          Cells[2,0]:='Direct '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[3,0]:='Direct '+flux_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

        end
        else begin
          Cells[2,0]:='Diffuse '+irradiance_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+j,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          if wantPhotonflux=1 then
            Cells[3,0]:='Diffuse '+flux_title;
        end;
    end;

    if coordinateSystem=2 then
    begin
      Cells[0,0]:='Inclination-Angle_(TTA)';

      if (DirectIrradiance) and (DiffuseIrradiance) then
      begin
        Cells[1,0]:='Direct '+irradiance_title;
        Cells[2,0]:='Diffuse '+irradiance_title;
        if wantPhotonData then
        begin
          for j := 0 to numberOfWavelengths-1 do
            Cells[3+2*wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
          for j := 0 to numberOfWavelengths-1 do
            Cells[3+2*wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
        if wantPhotonflux=1 then
        begin
          Cells[3,0]:='Direct '+flux_title;
          Cells[4,0]:='Diffuse '+flux_title;
        end;
      end else
        if DirectIrradiance then
        begin
          Cells[1,0]:='Direct '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[2,0]:='Direct '+flux_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+j,0]:='Direct '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';

        end
        else begin
          Cells[1,0]:='Diffuse '+irradiance_title;
          if wantPhotonflux=1 then
            Cells[2,0]:='Diffuse '+flux_title;
          if wantPhotonData then
            for j := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+j+numberOfWavelengths,0]:='Diffuse '+FloatToStr(startWavelength + j*intervalWavelength)+'nm [#Photons/(cm^2*nm)]';
        end;
    end;

    // data

    if (coordinateSystem=0) or (coordinateSystem=1) or (coordinateSystem=3) then
      for j := 1 to RowCount do
      begin
        Cells[0,j]:=FloatToStr(PostProcessingMatrix[1,j]); // AZI bzw. YAN
        Cells[1,j]:=FloatToStr(PostProcessingMatrix[2,j]); // ALT bzw. ZAN

        if DirectIrradiance and DiffuseIrradiance then
        begin
          Cells[2,j]:=FloatToStr(PostProcessingMatrix[3,j]); // Direct irradiance
          Cells[3,j]:=FloatToStr(PostProcessingMatrix[4,j]); // Diffuse irradiance
          if wantPhotonflux=1 then
          begin
            Cells[4,j]:=FloatToStr(PostProcessingMatrix[5,j]); // Direct photon flux
            Cells[5,j]:=FloatToStr(PostProcessingMatrix[6,j]); // Diffuse photon flux
          end;
          if wantPhotonData then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+counter,j]:=FloatToStr(PostProcessingMatrix[7+counter,j]);
            for counter := 0 to numberOfWavelengths-1 do
              Cells[4+2*wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(PostProcessingMatrix[7+numberOfWavelengths+counter,j]);
          end;
          
        end else
        if DirectIrradiance then
        begin
          Cells[2,j]:=FloatToStr(PostProcessingMatrix[3,j]); // Direct irradiance
          if wantPhotonflux=1 then
            Cells[3,j]:=FloatToStr(PostProcessingMatrix[5,j]); // Direct photon flux
          if wantPhotonData then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+counter,j]:=FloatToStr(PostProcessingMatrix[7+counter,j]);
        end
        else begin
          Cells[2,j]:=FloatToStr(PostProcessingMatrix[4,j]); // Diffuse irradiance
          if wantPhotonflux=1 then
            Cells[3,j]:=FloatToStr(PostProcessingMatrix[6,j]); // Diffuse photon flux
          if wantPhotonData then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(PostProcessingMatrix[7+numberOfWavelengths+counter,j]);
        end;
      end;

    if (coordinateSystem=2) then
      for j := 1 to RowCount do
      begin
        Cells[0,j]:=FloatToStr(PostProcessingMatrix[1,j]); // Matrix export

        if DirectIrradiance and DiffuseIrradiance then
        begin
          Cells[1,j]:=FloatToStr(PostProcessingMatrix[3,j]); // Direct irradiance
          Cells[2,j]:=FloatToStr(PostProcessingMatrix[4,j]); // Diffuse irradiance
          if wantPhotonflux=1 then
          begin
            Cells[3,j]:=FloatToStr(PostProcessingMatrix[5,j]); // Direct photon flux
            Cells[4,j]:=FloatToStr(PostProcessingMatrix[6,j]); // Diffuse photon flux
          end;
          if wantPhotonData then
          begin
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+counter,j]:=FloatToStr(PostProcessingMatrix[7+counter,j]);
            for counter := 0 to numberOfWavelengths-1 do
              Cells[3+2*wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(PostProcessingMatrix[7+numberOfWavelengths+counter,j]);
          end;
        end else
        if DirectIrradiance then
        begin
          Cells[1,j]:=FloatToStr(PostProcessingMatrix[3,j]); // Direct irradiance
          if wantPhotonflux=1 then
            Cells[2,j]:=FloatToStr(PostProcessingMatrix[5,j]); // Direct photon flux
          if wantPhotonData then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+counter,j]:=FloatToStr(PostProcessingMatrix[7+counter,j]);
        end
        else begin
          Cells[1,j]:=FloatToStr(PostProcessingMatrix[4,j]); // Diffuse irradiance
          if wantPhotonflux=1 then
            Cells[2,j]:=FloatToStr(PostProcessingMatrix[6,j]); // Diffuse photon flux
          if wantPhotonData then
            for counter := 0 to numberOfWavelengths-1 do
              Cells[2+wantPhotonflux+counter+numberOfWavelengths,j]:=FloatToStr(PostProcessingMatrix[7+numberOfWavelengths+counter,j]);
        end;
      end;
  end;  // with OutputGrid do begin
end;




function calculateYearlyIrradianceBinMatrix (radiationType : TRadiationType): double; STDCALL;
var
  radiationSum : double;
  i : integer;
begin
  if isBinMatrixInitialized then
  begin
    radiationSum := 0;

    for i := 1 to 129600 do
    begin
      if radiationType=[direct] then
        radiationSum := radiationSum + BinMatrix[3,i] // Direct
      else
        radiationSum := radiationSum + BinMatrix[4,i]; // Diffuse
    end;

    result := radiationSum;
  end else
    result := -1;
end;

function calculateYearlyIrradiancePostProcessingMatrix (radiationType : TRadiationType): double; STDCALL;
var
  radiationSum : double;
  i : integer;
begin
  if isPostProcessingMatrixInitialized then
  begin
    radiationSum := 0;

    for i := 1 to 129600 do
    begin
      if radiationType=[direct] then
        radiationSum := radiationSum + PostProcessingMatrix[3,i] // Direct
      else
        radiationSum := radiationSum + PostProcessingMatrix[4,i]; // Diffuse
    end;

    result := radiationSum;
  end else
    result := -1;
end;


function normalizePostProcessingMatrix(directTargetIrradiance,diffuseTargetIrradiance : double; wantPhotonData:boolean; numberOfWavelengths : integer):boolean; STDCALL;
var
  directYearlyIrradiance, diffuseYearlyIrradiance,
  globalYearlyIrradiance : double;
  scalingFactor_direct,scalingFactor_diffuse : double;
  i,counter : integer;
begin

  directYearlyIrradiance := calculateYearlyIrradianceBinMatrix([direct])/(1000*60*60);
  diffuseYearlyIrradiance := calculateYearlyIrradianceBinMatrix([diffuse])/(1000*60*60);

  if ((directYearlyIrradiance<0) or (diffuseYearlyIrradiance<0)) then
    result := false
  else begin
    globalYearlyIrradiance := directYearlyIrradiance + diffuseYearlyIrradiance;
    scalingFactor_direct := directTargetIrradiance / directYearlyIrradiance;
    scalingFactor_diffuse := diffuseTargetIrradiance / diffuseYearlyIrradiance;


    for i := 1 to 129600 do
    begin
      PostProcessingMatrix[1,i] := BinMatrix[1,i];
      PostProcessingMatrix[2,i] := BinMatrix[2,i];
      PostProcessingMatrix[3,i] := BinMatrix[3,i] * scalingFactor_direct;
      PostProcessingMatrix[4,i] := BinMatrix[4,i] * scalingFactor_diffuse;

      PostProcessingMatrix[5,i] := BinMatrix[5,i] * scalingFactor_direct;
      PostProcessingMatrix[6,i] := BinMatrix[6,i] * scalingFactor_diffuse;

      for counter := 0 to numberOfWavelengths - 1 do
      begin
        PostProcessingMatrix[7 + counter,i]:=BinMatrix[7 + counter,i] * scalingFactor_direct;
        PostProcessingMatrix[7 + numberOfWavelengths + counter,i]:=BinMatrix[7 + numberOfWavelengths + counter,i] * scalingFactor_diffuse;
      end;
    end;
    result := true;
  end;

end;


exports
  resetInitialization,
  InitializeBinning,
  InitializePostProcessing,
  BinValues,
  BinningOutput,
  BinningOutputVar,
  PostProcessingOutput,
  calculateYearlyIrradianceBinMatrix,
  calculateYearlyIrradiancePostProcessingMatrix,
  normalizePostProcessingMatrix;

begin
end.