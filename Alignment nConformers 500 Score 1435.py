import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 12, 41623])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwxOfYdVCWJhbGxTY2FsZXEDSwxHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLDEc/8AAAAAAAAH2HVQVjb2xvcnEFSwxLAH1xBihLAV1xB0sBYUsCXXEISwJhSwNdcQlLA2FLBF1xCksEYUsFXXELSwVhSwZdcQxLBmFLB11xDUsHYUsIXXEOSwhhSwldcQ9LCWFLCl1xEEsKYUsLXXERSwthdYdVCnJpYmJvblR5cGVxEksMSwB9h1UKc3RpY2tTY2FsZXETSwxHP/AAAAAAAAB9h1UMbW1DSUZIZWFkZXJzcRRdcRUoTk5OTk5OTk5OTk5OZVUMYXJvbWF0aWNNb2RlcRZLDEsBfYdVCnZkd0RlbnNpdHlxF0sMR0AUAAAAAAAAfYdVBmhpZGRlbnEYSwyJfYdVDWFyb21hdGljQ29sb3JxGUsMTn2HVQ9yaWJib25TbW9vdGhpbmdxGksMSwB9h1UJYXV0b2NoYWlucRtLDIh9h1UKcGRiVmVyc2lvbnEcSwxLAH2HVQhvcHRpb25hbHEdfXEeVQhvcGVuZWRBc3EfiIlLDChVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzRrN2FfbGlnYW5kLnNkZnEgVQtNREwgTU9ML1NERnEhTol0cSJ9cSMoKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzJhbTlfbGlnYW5kLnNkZnEkaCFOiXRxJV1xJksBYShVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzF4cTNfbGlnYW5kLnNkZnEnaCFOiXRxKF1xKUsKYShVKi9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy9SXzFlM2tfbGlnYW5kLnNkZnEqaCFOiXRxK11xLEsJYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8yaHZjX2xpZ2FuZC5zZGZxLWghTol0cS5dcS9LBmEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvMmF4OV9saWdhbmQuc2RmcTBoIU6JdHExXXEySwVhKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzJhbWJfbGlnYW5kLnNkZnEzaCFOiXRxNF1xNUsDYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8yYXg2X2xpZ2FuZC5zZGZxNmghTol0cTddcThLBGEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvM3Jsal9saWdhbmQuc2RmcTloIU6JdHE6XXE7SwdhKFUoL2hvbWUvaW5zaWxpY2hlbS9MaWdhbmRzLzFlM2dfbGlnYW5kLnNkZnE8aCFOiXRxPV1xPksAYShVKC9ob21lL2luc2lsaWNoZW0vTGlnYW5kcy8zcmxsX2xpZ2FuZC5zZGZxP2ghTol0cUBdcUFLCGEoVSgvaG9tZS9pbnNpbGljaGVtL0xpZ2FuZHMvMmFtYV9saWdhbmQuc2RmcUJoIU6JdHFDXXFESwJhdYeHc1UPbG93ZXJDYXNlQ2hhaW5zcUVLDIl9h1UJbGluZVdpZHRocUZLDEc/8AAAAAAAAH2HVQ9yZXNpZHVlTGFiZWxQb3NxR0sMSwB9h1UEbmFtZXFISwxYDwAAADJIVkNfTEdEX0FfMjIyNn1xSShYDwAAADFFM0tfUjE4X0FfMTAwMF1xSksJYVgPAAAAMUUzR19SMThfQV8xMDAwXXFLSwBhWAwAAAAzUkxMX1JMTF9BXzFdcUxLCGFYDwAAADFYUTNfUjE4X0FfMTAwMV1xTUsKYVgMAAAAMkFYOV9CSE1fQV8xXXFOSwVhWAwAAAAzUkxKX1JMSl9BXzFdcU9LB2FYDwAAADJBTUFfREhUX0FfMTAwMV1xUEsCYVgPAAAAMkFNQl8xN0hfQV8xMDAxXXFRSwNhWA8AAAAyQU05X1RFU19BXzEwMDBdcVJLAWFYDwAAADRLN0FfREhUX0FfMTAwMV1xU0sLYVgMAAAAMkFYNl9IRlRfQV8xXXFUSwRhdYdVD2Fyb21hdGljRGlzcGxheXFVSwyJfYdVD3JpYmJvblN0aWZmbmVzc3FWSwxHP+mZmZmZmZp9h1UKcGRiSGVhZGVyc3FXXXFYKH1xWX1xWn1xW31xXH1xXX1xXn1xX31xYH1xYX1xYn1xY31xZGVVA2lkc3FlSwxLCUsAhn1xZihLAEsAhl1xZ0sAYUsHSwCGXXFoSwdhSwNLAIZdcWlLA2FLCEsAhl1xaksIYUsGSwCGXXFrSwZhSwtLAIZdcWxLC2FLAksAhl1xbUsCYUsFSwCGXXFuSwVhSwpLAIZdcW9LCmFLAUsAhl1xcEsBYUsESwCGXXFxSwRhdYdVDnN1cmZhY2VPcGFjaXR5cXJLDEe/8AAAAAAAAH2HVRBhcm9tYXRpY0xpbmVUeXBlcXNLDEsCfYdVFHJpYmJvbkhpZGVzTWFpbmNoYWlucXRLDIh9h1UHZGlzcGxheXF1SwyIfYd1Lg=='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksMVQEgfYdVC2ZpbGxEaXNwbGF5cQNLDIl9h1UEbmFtZXEESwxYAwAAAFVOS32HVQVjaGFpbnEFSwxYAQAAACB9h1UOcmliYm9uRHJhd01vZGVxBksMSwJ9h1UCc3NxB0sMiYmGfYdVCG1vbGVjdWxlcQhLDEsAfXEJKEsBTl1xCksBSwGGcQthhksCTl1xDEsCSwGGcQ1hhksDTl1xDksDSwGGcQ9hhksETl1xEEsESwGGcRFhhksFTl1xEksFSwGGcRNhhksGTl1xFEsGSwGGcRVhhksHTl1xFksHSwGGcRdhhksITl1xGEsISwGGcRlhhksJTl1xGksJSwGGcRthhksKTl1xHEsKSwGGcR1hhksLTl1xHksLSwGGcR9hhnWHVQtyaWJib25Db2xvcnEgSwxOfYdVBWxhYmVscSFLDFgAAAAAfYdVCmxhYmVsQ29sb3JxIksMTn2HVQhmaWxsTW9kZXEjSwxLAX2HVQVpc0hldHEkSwyJfYdVC2xhYmVsT2Zmc2V0cSVLDE59h1UIcG9zaXRpb25xJl1xJyhLAUsBhnEoSwFLAYZxKUsBSwGGcSpLAUsBhnErSwFLAYZxLEsBSwGGcS1LAUsBhnEuSwFLAYZxL0sBSwGGcTBLAUsBhnExSwFLAYZxMksBSwGGcTNlVQ1yaWJib25EaXNwbGF5cTRLDIl9h1UIb3B0aW9uYWxxNX1VBHNzSWRxNksMSv////99h3Uu'))
	atomInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQdyZXNpZHVlcQJNFAFLFH1xAyhLDE5dcQRLAEsVhnEFYYZLDU5dcQZLFUsVhnEHYYZLDk5dcQhLKksVhnEJYYZLD05dcQpLP0sXhnELYYZLEE5dcQxLVksUhnENYYZLEU5dcQ5LaksVhnEPYYZLEk5dcRBLf0sahnERYYZLE05dcRJLmUschnETYYZLFU5dcRRL1UsVhnEVYYZLFk5dcRZL6ksVhnEXYYZLF05dcRhL/0sVhnEZYYZ1h1UIdmR3Q29sb3JxGk0UAU59h1UEbmFtZXEbTRQBWAMAAABDMTF9cRwoWAMAAABDMTldcR0oSxJLKUs+S1FLs0vQS+dL/E0TAWVYAwAAAEMxOF1xHihLEUsoSz1LUEuyS89L5kv7TRIBZVgDAAAAQzEzXXEfKEsMSyJLN0tLS5NLrUvIS+FL9k0MAWVYAwAAAEMxMl1xIChLC0shSzZLSkuPS6tLx0vgS/VNCwFlWAMAAABDMTBdcSEoSwlLH0s0S0hLZ0t6S4xLqEvFS95L800JAWVYAwAAAEMxN11xIihLEEsmSztLT0uxS85L5Uv6TRABZVgCAAAATzJdcSMoSxRLJ0s8S1VLX0twS6lLuUvpS/5NEQFlWAIAAABPMV1xJChLE0sYSy1LVEteS21LmEunS7hL6Ev9TQIBZVgDAAAAQzE0XXElKEsNSyNLOEtMS5dLrkvLS+JL900NAWVYAgAAAE80XXEmKEtpS3llWAIAAABPM11xJyhLZUt4S6xLymVYAwAAAEJyMV1xKEtqYVgDAAAAQzE2XXEpKEsPSyVLOktOS7BLzUvkS/lNDwFlWAIAAABDOV1xKihLCEseSzNLR0tmS3ZLi0umS8RL3UvyTQgBZVgCAAAAQzhdcSsoSwdLHUsyS0ZLZEt1S4pLo0vDS9xL8U0HAWVYAwAAAEMxNV1xLChLDkskSzlLTUuvS8xL40v4TQ4BZVgCAAAAQzNdcS0oSwJLF0ssS0FLW0tuS4JLnUu+S9dL7E0BAWVYAgAAAEMyXXEuKEsBSxZLK0tAS1pLbEuBS5tLvUvWS+tNAAFlWAIAAABDMV1xLyhLAEsVSypLP0tXS2tLgEuZS7VL1UvqS/9lWAIAAABDN11xMChLBkscSzFLRUtiS3RLiUuiS8JL20vwTQYBZVgCAAAAQzZdcTEoSwVLG0swS0RLYUtzS4hLoUvBS9pL700FAWVYAgAAAEM1XXEyKEsESxpLL0tDS2BLckuHS6BLwEvZS+5NBAFlWAIAAABDNF1xMyhLA0sZSy5LQktcS29LhkufS79L2EvtTQMBZVgDAAAAQzIyXXE0S9NhWAMAAABDMjNdcTVL1GFYAwAAAEMyMF1xNihLUkvRZVgDAAAAQzIxXXE3KEtTS9JlWAIAAABOMV1xOChLXUtxS39LpEu2ZVgCAAAATjJdcTkoS2NLd0uOS6VLt2VYAgAAAE4zXXE6KEu0S8llWAIAAABGMV1xOyhLVkt8S4NLmku6ZVgCAAAARjJdcTwoS1hLfUuES5xLu2VYAgAAAEYzXXE9KEtZS35LhUueS7xlWAIAAABGNF1xPkuQYVgCAAAARjVdcT9LkWFYAgAAAEY2XXFAS5JhWAIAAABGN11xQUuUYVgCAAAARjhdcUJLlWFYAgAAAEY5XXFDS5ZhdYdVA3Zkd3FETRQBiX2HVQ5zdXJmYWNlRGlzcGxheXFFTRQBiX2HVQVjb2xvcnFGTRQBTn1xRyhLD11xSEtqYUsMXXFJKEsTSxRLGEsnSy1LPEtUS1VLXktfS2VLaUttS3BLeEt5S5hLp0upS6xLuEu5S8pL6EvpS/1L/k0CAU0RAWVLDV1xSihLVktYS1lLfEt9S35Lg0uES4VLkEuRS5JLlEuVS5ZLmkucS55Luku7S7xlSw5dcUsoS11LY0txS3dLf0uOS6RLpUu0S7ZLt0vJZXWHVQlpZGF0bVR5cGVxTE0UAYl9h1UGYWx0TG9jcU1NFAFVAH2HVQVsYWJlbHFOTRQBWAAAAAB9h1UOc3VyZmFjZU9wYWNpdHlxT00UAUe/8AAAAAAAAH2HVQdlbGVtZW50cVBNFAFLBn1xUShLCF1xUihLE0sUSxhLJ0stSzxLVEtVS15LX0tlS2lLbUtwS3hLeUuYS6dLqUusS7hLuUvKS+hL6Uv9S/5NAgFNEQFlSwldcVMoS1ZLWEtZS3xLfUt+S4NLhEuFS5BLkUuSS5RLlUuWS5pLnEueS7pLu0u8ZUsjXXFUS2phSwddcVUoS11LY0txS3dLf0uOS6RLpUu0S7ZLt0vJZXWHVQpsYWJlbENvbG9ycVZNFAFOfYdVDHN1cmZhY2VDb2xvcnFXTRQBTn2HVQ9zdXJmYWNlQ2F0ZWdvcnlxWE0UAVgEAAAAbWFpbn2HVQZyYWRpdXNxWU0UAUc//hR64AAAAH1xWihHP/nCj2AAAABdcVsoSwJLBEsISwlLF0saSyxLQ0tGS0dLSUtaS1xLYktkS29Lckt2S3pLgEuCS4ZLiEuLS51Ln0ugS6ZLrUuwS8tLzEvNS85Lz0vQS9FL0kvXS9lL3UveS+xL7kvyS/NNAQFlRz/2uFHgAAAAXXFcKEsTSxhLLUtUS15LX0tlS3BLeEt5S5hLp0u4S+hL/U0CAWVHP/dcKQAAAABdcV0oSxRLJ0s8S1VLaUttS6lLrEu5S8pL6Uv+TREBZUdAA1wpAAAAAF1xXktqYUc/9HrhQAAAAF1xXyhLVktYS1lLfEt9S35Lg0uES4VLkEuRS5JLlEuVS5ZLmkucS55Luku7S7xlRz/6PXCgAAAAXXFgKEtdS2NLcUt3S39LjkukS6VLtEu2S7dLyWVHP/wo9cAAAABdcWEoSwNLCksLSxlLQEtES0pLW0tgS2FLc0t0S3VLgUuHS4lLikuZS5tLoUuuS69LsUuyS79LwEvBS8JLw0vES8VLxkvHS9hL30vgS+1L9Ev1ZXWHVQpjb29yZEluZGV4cWJdcWMoSwBLFYZxZEsASxWGcWVLAEsVhnFmSwBLF4ZxZ0sASxSGcWhLAEsVhnFpSwBLGoZxaksASxyGcWtLAEsghnFsSwBLFYZxbUsASxWGcW5LAEsVhnFvZVULbGFiZWxPZmZzZXRxcE0UAU59h1USbWluaW11bUxhYmVsUmFkaXVzcXFNFAFHAAAAAAAAAAB9h1UIZHJhd01vZGVxck0UAUsCfYdVCG9wdGlvbmFscXN9cXQoVQxzZXJpYWxOdW1iZXJxdYiJTRQBSwF9cXYoSwJdcXcoSwFLFksrS0BLV0trS4BLmku2S9ZL600AAWVLA11xeChLAksXSyxLQUtYS2xLgUubS7dL10vsTQEBZUsEXXF5KEsDSxhLLUtCS1lLbUuCS5xLuEvYS+1NAgFlSwVdcXooSwRLGUsuS0NLWktuS4NLnUu5S9lL7k0DAWVLBl1xeyhLBUsaSy9LREtbS29LhEueS7pL2kvvTQQBZUsHXXF8KEsGSxtLMEtFS1xLcEuFS59Lu0vbS/BNBQFlSwhdcX0oSwdLHEsxS0ZLXUtxS4ZLoEu8S9xL8U0GAWVLCV1xfihLCEsdSzJLR0teS3JLh0uhS71L3UvyTQcBZUsKXXF/KEsJSx5LM0tIS19Lc0uIS6JLvkveS/NNCAFlSwtdcYAoSwpLH0s0S0lLYEt0S4lLo0u/S99L9E0JAWVLDF1xgShLC0sgSzVLSkthS3VLikukS8BL4Ev1TQoBZUsNXXGCKEsMSyFLNktLS2JLdkuLS6VLwUvhS/ZNCwFlSw5dcYMoSw1LIks3S0xLY0t3S4xLpkvCS+JL900MAWVLD11xhChLDksjSzhLTUtkS3hLjUunS8NL40v4TQ0BZUsQXXGFKEsPSyRLOUtOS2VLeUuOS6hLxEvkS/lNDgFlSxFdcYYoSxBLJUs6S09LZkt6S49LqUvFS+VL+k0PAWVLEl1xhyhLEUsmSztLUEtnS3tLkEuqS8ZL5kv7TRABZUsTXXGIKEsSSydLPEtRS2hLfEuRS6tLx0vnS/xNEQFlSxRdcYkoSxNLKEs9S1JLaUt9S5JLrEvIS+hL/U0SAWVLFV1xiihLFEspSz5LU0t+S5NLrUvJS+lL/k0TAWVLFl1xiyhLVEuUS65LymVLF11xjChLVUuVS69Ly2VLGF1xjShLlkuwS8xlSxldcY4oS5dLsUvNZUsaXXGPKEuYS7JLzmVLG11xkChLs0vPZUscXXGRKEu0S9BlSx1dcZJL0WFLHl1xk0vSYUsfXXGUS9NhSyBdcZVL1GF1h4dVB2JmYWN0b3JxloiJTRQBRwAAAAAAAAAAfYeHVQlvY2N1cGFuY3lxl4iJTRQBRz/wAAAAAAAAfYeHdVUHZGlzcGxheXGYTRQBiH2HdS4='))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECTS0BTn2HVQVhdG9tc3EDXXEEKF1xBShLGEsZZV1xBihLGEshZV1xByhLGUsaZV1xCChLGksbZV1xCShLGksrZV1xCihLG0scZV1xCyhLHEsdZV1xDChLHEshZV1xDShLHUseZV1xDihLHksfZV1xDyhLH0sgZV1xEChLH0slZV1xEShLIEshZV1xEihLIEsiZV1xEyhLIksjZV1xFChLI0skZV1xFShLJEslZV1xFihLJEsoZV1xFyhLJEspZV1xGChLJUsmZV1xGShLJksnZV1xGihLJ0soZV1xGyhLKEsqZV1xHChLKEssZV1xHShLLUsuZV1xHihLLUs3ZV1xHyhLLksvZV1xIChLL0swZV1xIShLL0sxZV1xIihLMUsyZV1xIyhLMkszZV1xJChLMks3ZV1xJShLM0s0ZV1xJihLNEs1ZV1xJyhLNUs2ZV1xKChLNUs7ZV1xKShLNks3ZV1xKihLNks4ZV1xKyhLN0tBZV1xLChLOEs5ZV1xLShLOUs6ZV1xLihLOks7ZV1xLyhLOks+ZV1xMChLOktAZV1xMShLO0s8ZV1xMihLPEs9ZV1xMyhLPUs+ZV1xNChLPks/ZV1xNShLQktDZV1xNihLQktMZV1xNyhLQ0tEZV1xOChLREtFZV1xOShLREtGZV1xOihLRktHZV1xOyhLR0tIZV1xPChLR0tMZV1xPShLSEtJZV1xPihLSUtKZV1xPyhLSktLZV1xQChLSktQZV1xQShLS0tMZV1xQihLS0tNZV1xQyhLTEtWZV1xRChLTUtOZV1xRShLTktPZV1xRihLT0tQZV1xRyhLT0tTZV1xSChLT0tVZV1xSShLUEtRZV1xSihLUUtSZV1xSyhLUktTZV1xTChLU0tUZV1xTShLV0tdZV1xTihLV0teZV1xTyhLWEteZV1xUChLWEtfZV1xUShLWUtaZV1xUihLWUtgZV1xUyhLWktfZV1xVChLW0tcZV1xVShLW0tgZV1xVihLW0thZV1xVyhLXEtiZV1xWChLXUthZV1xWShLXktsZV1xWihLX0thZV1xWyhLYEtkZV1xXChLYktjZV1xXShLY0tkZV1xXihLY0tnZV1xXyhLY0toZV1xYChLZEtlZV1xYShLZUtmZV1xYihLZktnZV1xYyhLZ0trZV1xZChLZ0ttZV1xZShLaEtpZV1xZihLaktrZV1xZyhLbktvZV1xaChLb0twZV1xaShLb0txZV1xaihLb0tyZV1xayhLcktzZV1xbChLckt0ZV1xbShLc0t6ZV1xbihLdEt1ZV1xbyhLdEt4ZV1xcChLdUt2ZV1xcShLdUt3ZV1xcihLeEt5ZV1xcyhLeUt6ZV1xdChLekt7ZV1xdShLe0t8ZV1xdihLfEt9ZV1xdyhLfEt+ZV1xeChLfkt/ZV1xeShLfkuAZV1xeihLfkuBZV1xeyhLgkuDZV1xfChLg0uEZV1xfShLhEuFZV1xfihLhEuGZV1xfyhLhEuHZV1xgChLh0uIZV1xgShLh0uJZV1xgihLiUuKZV1xgyhLikuLZV1xhChLikuMZV1xhShLi0uSZV1xhihLjEuNZV1xhyhLjUuOZV1xiChLjkuPZV1xiShLjkuSZV1xiihLj0uQZV1xiyhLj0uRZV1xjChLkkuTZV1xjShLk0uUZV1xjihLk0uVZV1xjyhLk0uWZV1xkChLl0uYZV1xkShLl0ujZV1xkihLmEuwZV1xkyhLmEuZZV1xlChLmUuaZV1xlShLmkueZV1xlihLmkukZV1xlyhLm0ukZV1xmChLnEukZV1xmShLnUukZV1xmihLnkujZV1xmyhLnkufZV1xnChLn0ugZV1xnShLoEuhZV1xnihLoEumZV1xnyhLoUuiZV1xoChLokujZV1xoShLpUumZV1xoihLpUurZV1xoyhLpkunZV1xpChLp0uvZV1xpShLqEurZV1xpihLqUurZV1xpyhLqkurZV1xqChLrEuvZV1xqShLrUuvZV1xqihLrkuvZV1xqyhLsUuzZV1xrChLsUu4ZV1xrShLsku6ZV1xrihLs0u1ZV1xryhLtEu6ZV1xsChLtUu3ZV1xsShLtUu7ZV1xsihLtku6ZV1xsyhLt0u5ZV1xtChLt0u6ZV1xtShLuEu5ZV1xtihLuEu9ZV1xtyhLu0u8ZV1xuChLvUu+ZV1xuShLvku/ZV1xuihLvkvAZV1xuyhLwEvBZV1xvChLwEvCZV1xvShLwEvDZV1xvihLw0vEZV1xvyhLxEvFZV1xwChLxUvGZV1xwShLxUvKZV1xwihLxkvHZV1xwyhLx0vIZV1xxChLyEvJZV1xxShLyEvLZV1xxihLyUvKZV1xxyhLy0vMZV1xyChLzUvrZV1x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'))
	crdInfo = cPickle.loads(base64.b64decode('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'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {u'': ((0.921569, 0, 0.14902), 1, u''), u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'gold': ((1, 0.843137, 0), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), u'Pr': ((0.85098, 1, 0.780392), 1, u'default'),
u'deep pink': ((1, 0.0784314, 0.576471), 1, u'default'), u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), u'Pu': ((0, 0.419608, 1), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'purple': ((0.627451, 0.12549, 0.941176), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'),
u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), u'sky blue': ((0.529412, 0.807843, 0.921569), 1, u'default'), u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'Mt': ((0.921569, 0, 0.14902), 1, u'default'), u'plum': ((0.866667, 0.627451, 0.866667), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'), u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'dodger blue': ((0.117647, 0.564706, 1), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'),
u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Tm': ((0, 0.831373, 0.321569), 1, u'default'), u'light green': ((0.564706, 0.933333, 0.564706), 1, u'default'),
u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), u'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'salmon': ((0.980392, 0.501961, 0.447059), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'),
u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'), u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default'), u'light gray': ((0.827451, 0.827451, 0.827451), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30), u'': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 16, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (18, (u'deep pink', (1, 0.0784314, 0.576471, 1)), {(u'green', (0, 1, 0, 1)): [17], (u'Br', (0.65098, 0.160784, 0.160784, 1)): [15], (u'light green', (0.564706, 0.933333, 0.564706, 1)): [3], (u'dodger blue', (0.117647, 0.564706, 1, 1)): [8], (u'F', (0.564706, 0.878431, 0.313725, 1)): [13], (u'', (0.358033, 0.260402, 0.804281, 1)): [10], (u'N', (0.188235, 0.313725, 0.972549, 1)): [14], (u'', (0.123026, 0.337066, 0.083208, 1)): [11], (u'purple', (0.627451, 0.12549, 0.941176, 1)): [9], (u'gold', (1, 0.843137, 0, 1)): [7], (u'sky blue', (0.529412, 0.807843, 0.921569, 1)): [1], (u'tan', (0.823529, 0.705882, 0.54902, 1)): [0], (u'O', (1, 0.0509804, 0.0509804, 1)): [12], (u'plum', (0.866667, 0.627451, 0.866667, 1)): [2], (u'light gray', (0.827451, 0.827451, 0.827451, 1)): [5], (u'salmon', (0.980392, 0.501961, 0.447059, 1)): [4], (u'yellow', (1, 1, 0, 1)): [16]})
	viewerInfo = {'cameraAttrs': {'center': (0.10849997615814, 30.985500019073, 4.3005), 'fieldOfView': 25.350031531442, 'nearFar': (16.079583727155, -12.694803936947), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 4.3005}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 10.357999968211, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 1, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 17, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(12, 'Chimera default', 'rounded', u'unknown'), (13, 'Chimera default', 'rounded', u'unknown'), (14, 'Chimera default', 'rounded', u'unknown'), (15, 'Chimera default', 'rounded', u'unknown'), (16, 'Chimera default', 'rounded', u'unknown'), (17, 'Chimera default', 'rounded', u'unknown'), (18, 'Chimera default', 'rounded', u'unknown'), (19, 'Chimera default', 'rounded', u'unknown'), (20, 'Chimera default', 'rounded', u'unknown'), (21, 'Chimera default', 'rounded', u'unknown'), (22, 'Chimera default', 'rounded', u'unknown'), (23, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.3574067443365933, 0.6604015517481455, -0.6604015517481456), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.3574067443365933, 0.6604015517481455, 0.6604015517481456), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.2505628070857316, 0.2505628070857316, 0.9351131265310294), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (512, 384)
	xformMap = {0: (((0, 0, 1), 0), (0, 0, 0), True), 1: (((0, 0, 1), 0), (0, 0, 0), True), 2: (((0, 0, 1), 0), (0, 0, 0), True), 3: (((0, 0, 1), 0), (0, 0, 0), True), 4: (((0, 0, 1), 0), (0, 0, 0), True), 5: (((0, 0, 1), 0), (0, 0, 0), True), 6: (((0, 0, 1), 0), (0, 0, 0), True), 7: (((0, 0, 1), 0), (0, 0, 0), True), 8: (((0, 0, 1), 0), (0, 0, 0), True), 9: (((0, 0, 1), 0), (0, 0, 0), True), 10: (((0, 0, 1), 0), (0, 0, 0), True), 11: (((0, 0, 1), 0), (0, 0, 0), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True, 8: True, 9: True, 10: True, 11: True, 601: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

