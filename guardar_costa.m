m_proj('lambert','lon',[min(min(double(lon))) max(max(double(lon)))],...
        'lat',[min(min(double(lat))) max(max(double(lat)))]);
  % set up projection parameters
  
  % This command does not draw anything - it merely processes the 
  % high-resolution database using the current projection parameters 
  % to generate a smaller coastline file called "gumby"
%   m_gshhs('fb1');             % Full resolution national borders

  m_gshhs_h('save','gumby');
  
  % Now we can draw a few maps of the same area much more quickly
  
  figure(1);
  m_usercoast('gumby','patch','r');
  m_grid;
  
  figure(2);

  m_grid('tickdir','out','yaxisloc','left');
  
  etc.