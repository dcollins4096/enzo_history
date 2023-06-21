;  Converts binary data into Katrin's "cosmo" dataset format.  This assumes
;  that the original dataset has already been preprocessed by GetPartPos_Binary
;  Note that we do some swap_endian stuff down there - Katrin wants
;  the cosmo stuff in little-endian, so if this is being run on a
;  little-endian machine there is no need to do any swapping.
pro convdm 

  infilename = 'unigrid1024/DD0099/dmposinfo.dat'
  outfilename = 'unigrid_1024_z_0.cosmo'

  print, 'infile is ', infilename
  print, 'outfile is ', outfilename

  edgesize = 0LL
  totalsize = 0LL

  edgesize = 256
  totalsize = 16777216
  ;totalsize = 10000

  close,1
  openr,1,infilename

  close,2
  openw,2,outfilename

  ;thisfloatval = 0.0

  xpos = 0.0 
  ypos = 0.0 
  zpos = 0.0 
  xvel = 0.0 
  yvel = 0.0 
  zvel = 0.0 
  mass = 0.0
  ;index = 0L

  minx = 1.0e30
  miny = 1.0e30
  minz = 1.0e30

  maxx = 0.0
  maxy = 0.0
  maxz = 0.0

  minxvel = 1.0e30
  minyvel = 1.0e30
  minzvel = 1.0e30

  maxxvel = -1.0e30 
  maxyvel = -1.0e30 
  maxzvel = -1.0e30

  omega_cdm = 0.27
  omega_b = 0.044
  omega_lambda = 0.686
  hubble = 0.71
  sigma_8 = 0.84
  n_s = 0.99
  z_start = 50.0
  z_now = 0.0
  comovingboxsize = 63.994   ; in Mpc/h

  lengthconvert = comovingboxsize / hubble
  rhozero = (omega_cdm + omega_b) * 1.8788e-29 * hubble * hubble
  timeconvert = 1.0 / (rhozero * 4.0 * 3.14159 * 6.67e-8 * (1.0+z_start)^3.0)^0.5
  velconvert = (lengthconvert*3.0824e24/timeconvert)*(1.0+z_now)/(1.0+z_start)/1.0e5;

  print, 'length conversion is ', lengthconvert
  print, 'time conversion is ', timeconvert
  print, 'velocity conversion is ', velconvert

  thisindex=0LL

  thisindex = 0

  print, 'looping over ', totalsize, ' values'

  while (thisindex lt totalsize) do begin

      ; read in floating point values: xpos, xvel, ypos, yvel,
      ; zpos, zvel

      ; ------- x position
      readu, 1, xpos
      xpos = xpos * lengthconvert

      if(xpos GT maxx) then maxx = xpos
      if(xpos LT minx) then minx = xpos

      xpos = swap_endian(xpos)


      ; ------- x velocity
      readu, 1, xvel
      xvel = xvel * velconvert

      if(xvel GT maxxvel) then maxxvel = xvel
      if(xvel LT minxvel) then minxvel = xvel

      xvel = swap_endian(xvel)

      ; ------- y position
      readu, 1, ypos
      ypos = ypos * lengthconvert

      if(ypos GT maxy) then maxy = ypos
      if(ypos LT miny) then miny = ypos

      ypos = swap_endian(ypos)

      ; ------- y velocity
      readu, 1, yvel
      yvel = yvel * velconvert

      if(yvel GT maxyvel) then maxyvel = yvel
      if(yvel LT minyvel) then minyvel = yvel

      yvel = swap_endian(yvel)

      ; ------- z position
      readu, 1, zpos
      zpos = zpos * lengthconvert

      if(zpos GT maxz) then maxz = zpos
      if(zpos LT minz) then minz = zpos

      zpos = swap_endian(zpos)

      ; ------- z velocity
      readu, 1, zvel
      zvel = zvel * velconvert

      if(zvel GT maxzvel) then maxzvel = zvel
      if(zvel LT minzvel) then minzvel = zvel

      zvel = swap_endian(zvel)

      ; -------- mass and index number!
      mass = 1.0
      mass = swap_endian(mass)

      index = Long(thisindex)
      index = swap_endian(index)

      ;printf, 2, xpos, xvel, ypos, yvel, zpos, zvel, mass, index
      writeu, 2, xpos, xvel, ypos, yvel, zpos, zvel, mass, index

      thisindex = thisindex+1
      
  endwhile

  print, 'min, max x pos', minx, maxx
  print, 'min, max y pos', miny, maxy
  print, 'min, max z pos', minz, maxz


  print, 'min, max x vel', minxvel, maxxvel
  print, 'min, max y vel', minyvel, maxyvel
  print, 'min, max z vel', minzvel, maxzvel

  close, 1
  close, 2

end
