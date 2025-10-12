"""
Compose comparison panels for LIVE 2.0 cluster vs PubChem molecule.
"""
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path


def compose_panel(
    left_img_path: str,
    right_img_path: str,
    title_left: str,
    title_right: str,
    out_path: str,
    subtitle_left: str = "",
    subtitle_right: str = ""
):
    """
    Create a side-by-side comparison panel.
    
    Args:
        left_img_path: Path to left image (LIVE 2.0 cluster)
        right_img_path: Path to right image (PubChem molecule)
        title_left: Title for left image
        title_right: Title for right image
        out_path: Output path for composed panel
        subtitle_left: Optional subtitle for left image
        subtitle_right: Optional subtitle for right image
    """
    # Load images
    left = Image.open(left_img_path).convert("RGBA")
    right = Image.open(right_img_path).convert("RGBA")
    
    # Unify height to 640px
    target_height = 640
    
    def scale_to_height(img: Image.Image, target_h: int) -> Image.Image:
        """Scale image to target height while preserving aspect ratio."""
        ratio = target_h / img.height
        new_width = int(img.width * ratio)
        return img.resize((new_width, target_h), Image.LANCZOS)
    
    left = scale_to_height(left, target_height)
    right = scale_to_height(right, target_height)
    
    # Layout parameters
    gap = 40
    title_height = 60
    subtitle_height = 40 if (subtitle_left or subtitle_right) else 0
    header_height = title_height + subtitle_height
    footer_height = 40
    
    canvas_width = left.width + right.width + gap * 3
    canvas_height = header_height + target_height + footer_height
    
    # Create dark canvas
    canvas = Image.new("RGBA", (canvas_width, canvas_height), (22, 22, 22, 255))
    draw = ImageDraw.Draw(canvas)
    
    # Try to load a better font
    try:
        # Try common font locations
        font_large = ImageFont.truetype("arial.ttf", 24)
        font_small = ImageFont.truetype("arial.ttf", 16)
    except:
        try:
            font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24)
            font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16)
        except:
            # Fallback to default
            font_large = ImageFont.load_default()
            font_small = ImageFont.load_default()
    
    # Draw titles
    draw.text((gap, 16), title_left, (255, 255, 255), font=font_large)
    draw.text((gap + left.width + gap, 16), title_right, (255, 255, 255), font=font_large)
    
    # Draw subtitles if provided
    if subtitle_left or subtitle_right:
        y_subtitle = title_height
        if subtitle_left:
            draw.text((gap, y_subtitle), subtitle_left, (200, 200, 200), font=font_small)
        if subtitle_right:
            draw.text((gap + left.width + gap, y_subtitle), subtitle_right, (200, 200, 200), font=font_small)
    
    # Paste images
    y_offset = header_height
    canvas.paste(left, (gap, y_offset), left)
    canvas.paste(right, (gap + left.width + gap, y_offset), right)
    
    # Draw separator line
    line_x = gap + left.width + gap // 2
    draw.line(
        [(line_x, y_offset), (line_x, y_offset + target_height)],
        fill=(80, 80, 80),
        width=2
    )
    
    # Save result
    canvas = canvas.convert("RGB")  # Convert to RGB before saving as PNG
    canvas.save(out_path)
    print(f"✓ Composed panel saved to: {out_path}")


def compose_panel_with_metadata(
    left_img_path: str,
    right_img_path: str,
    cluster_metadata: dict,
    pubchem_data: dict,
    out_path: str
):
    """
    Create a comparison panel with automatic title/subtitle generation.
    
    Args:
        left_img_path: Path to LIVE 2.0 cluster image
        right_img_path: Path to PubChem molecule image
        cluster_metadata: Cluster metadata dict
        pubchem_data: PubChem result dict (from chem.pubchem_similar_top)
        out_path: Output path
    """
    # Left side (LIVE 2.0 cluster)
    title_left = "LIVE 2.0 Cluster"
    subtitle_left = f"Size: {cluster_metadata.get('size', '?')} | " \
                   f"Bonds: {cluster_metadata.get('bonds', '?')}"
    
    # Right side (PubChem)
    cid = pubchem_data.get("cid", "–")
    name = pubchem_data.get("name", "unknown")
    formula = pubchem_data.get("formula", "")
    
    # Truncate long names
    if len(name) > 40:
        name = name[:37] + "..."
    
    title_right = f"PubChem CID {cid}"
    subtitle_right = f"{name}"
    if formula:
        subtitle_right += f" | {formula}"
    
    compose_panel(
        left_img_path,
        right_img_path,
        title_left,
        title_right,
        out_path,
        subtitle_left,
        subtitle_right
    )

